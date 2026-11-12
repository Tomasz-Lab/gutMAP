#!/usr/bin/env python3
"""Run the pipeline in batches that fit a disk budget.

    batches.py plan --ids samples.csv --input-dir INPUT_DIR --staging-dir STAGING_DIR \\
                    --disk-budget 20T [--out next_batch.csv]
    batches.py run  --ids samples.csv --input-dir INPUT_DIR --staging-dir STAGING_DIR \\
                    --disk-budget 20T [--n-batches N]

Each sample reserves --reserve-factor x its FASTQ size, estimated from runinfo as
2 x bases + 120 B x spots. If runinfo has no bases for it, the sample reserves
--input-factor x the size of its dir in --input-dir. A batch takes the samples in the order of --ids until the
next sample does not fit in --disk-budget. plan prints the next batch and writes
nothing unless --out is given.

run creates <staging-dir>/batch_NNNN/ with batch.csv, runs the pipeline with that dir
as outdir, writes samples.tsv, copies the logs of failed tasks to failed_logs/ and
runs `nextflow clean -f` to remove the work dir of the run. It repeats until no sample
is eligible or after --n-batches batches. Start run from the launch dir of the
pipeline, one run at a time. A batch without samples.tsv (its run stopped) is
finished first. On SIGTERM or SIGINT, run passes the signal to nextflow, finishes
the batch and stops.

Sample states, from <staging-dir>/batch_NNNN/{batch.csv,samples.tsv}:
  running         in a batch that has no samples.tsv yet
  done, failed    the state in samples.tsv of its last batch
  retry limit     state retry after --max-batches batches
  not downloaded  not in --input-dir
  eligible        all other samples: new samples and samples to retry

In samples.tsv: failed = a task failed with an exit code that is not transient, or
the size gate failed; retry = a task failed with a transient exit code, or an
endpoint did not run; done = all endpoints ran and the size gate passed. If the run
did not end normally, failed becomes retry: a stopped task can end with any code.

Python 3.7 or later, standard library only.
"""
import argparse
import csv
import glob
import os
import re
import shutil
import signal
import subprocess
from collections import Counter, defaultdict

RUN_ACCESSION = re.compile(r"[SED]RR\d+")
BATCH_NAME = re.compile(r"batch_(\d+)")
FASTQ_BYTES_PER_SPOT = 120
SIZE_UNITS = {"": 1, "K": 10**3, "M": 10**6, "G": 10**9, "T": 10**12}
STATES = ["running", "done", "failed", "retry limit", "not downloaded", "eligible"]
OK_STATUS = {"COMPLETED", "CACHED"}
TRANSIENT_EXITS = {"-", "12", "130", "137", "140", "143", "247"}
ENDPOINTS = {"SYLPH_TAXONOMY", "CALCULATE_TPM_AND_ANNOTATE", "FILTER_AND_ANNOTATE"}
# A sylph profile with only its header has 218 bytes.
EMPTY_SYLPH_BYTES = 218
TASK_LOGS = [".command.sh", ".command.err", ".command.out", ".command.log", ".exitcode"]
DEFAULT_PIPELINE = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "mtg_pipeline.nf")


# --- inputs ---

def parse_size(text):
    """'20T' -> 20 * 10**12 bytes."""
    match = re.fullmatch(r"(\d+(?:\.\d+)?)\s*([KMGT]?)B?", text.strip().upper())
    if not match:
        raise argparse.ArgumentTypeError(f"not a size: {text}")
    return int(float(match.group(1)) * SIZE_UNITS[match.group(2)])


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def file_size(path):
    return os.path.getsize(path) if os.path.isfile(path) else 0


def dir_bytes(path):
    return sum(entry.stat().st_size for entry in os.scandir(path) if entry.is_file()) if os.path.isdir(path) else 0


def subdirs(path):
    return set(os.listdir(path)) if os.path.isdir(path) else set()


def read_csv(path, delimiter=","):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def read_ids(path):
    """Return the run accessions of the sra_id column in file order, without duplicates."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        if "sra_id" not in (reader.fieldnames or []):
            raise SystemExit(f"{path}: no sra_id column")
        ids = [(row["sra_id"] or "").strip() for row in reader]
    return list(dict.fromkeys(sra_id for sra_id in ids if RUN_ACCESSION.fullmatch(sra_id)))


def read_fastq_sizes(runinfo_csv):
    """Return {sra_id: estimated FASTQ bytes} for the runinfo rows with bases."""
    if not os.path.exists(runinfo_csv):
        return {}
    return {row["Run"]: 2 * number(row["bases"]) + FASTQ_BYTES_PER_SPOT * number(row["spots"])
            for row in read_csv(runinfo_csv) if number(row["bases"]) > 0}


# --- batches and sample states ---

def read_batches(staging_dir):
    """Return [(batch dir, sample ids, {sra_id: state} or None while it runs)] in batch order."""
    batches = []
    for name in sorted(subdirs(staging_dir)):
        batch_dir = os.path.join(staging_dir, name)
        if not BATCH_NAME.fullmatch(name) or not os.path.exists(os.path.join(batch_dir, "batch.csv")):
            continue
        ids = [row["sra_id"] for row in read_csv(os.path.join(batch_dir, "batch.csv"))]
        samples_tsv = os.path.join(batch_dir, "samples.tsv")
        states = None
        if os.path.exists(samples_tsv):
            states = {row["sra_id"]: row["state"] for row in read_csv(samples_tsv, delimiter="\t")}
        batches.append((batch_dir, ids, states))
    return batches


def next_batch_name(staging_dir):
    numbers = [0]
    for name in subdirs(staging_dir):
        match = BATCH_NAME.fullmatch(name)
        if match:
            numbers.append(int(match.group(1)))
    return f"batch_{max(numbers) + 1:04d}"


def sample_states(ids, batches, input_dir, max_batches):
    """Return {sra_id: one of STATES}."""
    downloaded = subdirs(input_dir)
    running, last_state, n_batches = set(), {}, Counter()
    for _, batch_ids, states in batches:
        n_batches.update(batch_ids)
        if states is None:
            running.update(batch_ids)
        else:
            last_state.update({sra_id: states.get(sra_id, "retry") for sra_id in batch_ids})

    result = {}
    for sra_id in ids:
        state = last_state.get(sra_id)
        if sra_id in running:
            result[sra_id] = "running"
        elif state in ("done", "failed"):
            result[sra_id] = state
        elif state == "retry" and n_batches[sra_id] >= max_batches:
            result[sra_id] = "retry limit"
        elif sra_id not in downloaded:
            result[sra_id] = "not downloaded"
        else:
            result[sra_id] = "eligible"
    return result


def choose_batch(eligible, reserve_of, budget):
    """Return ([(sra_id, reserved bytes)], too large, no size) for the next batch."""
    batch, total, too_large, no_size = [], 0, [], []
    for sra_id in eligible:
        reserved = reserve_of(sra_id)
        if reserved == 0:
            no_size.append(sra_id)
        elif reserved > budget:
            too_large.append(sra_id)
        elif total + reserved > budget:
            break
        else:
            batch.append((sra_id, reserved))
            total += reserved
    return batch, too_large, no_size


def next_batch(args):
    """Return (sample states, batches, next batch, too large, no size)."""
    ids = read_ids(args.ids)
    batches = read_batches(args.staging_dir)
    states = sample_states(ids, batches, args.input_dir, args.max_batches)
    eligible = [sra_id for sra_id in ids if states[sra_id] == "eligible"]
    fastq_sizes = read_fastq_sizes(args.runinfo)

    def reserve_of(sra_id):
        if sra_id in fastq_sizes:
            return int(args.reserve_factor * fastq_sizes[sra_id])
        return int(args.input_factor * dir_bytes(os.path.join(args.input_dir, sra_id)))

    batch, too_large, no_size = choose_batch(eligible, reserve_of, args.disk_budget)
    return states, batches, batch, too_large, no_size


def print_batch(name, batch, too_large, no_size, budget):
    reserved = sum(size for _, size in batch)
    print(f"{name}: {len(batch)} samples, reserve {reserved / 1e9:,.0f} GB of {budget / 1e9:,.0f} GB", flush=True)
    skipped_samples = (("reserve larger than the budget", too_large), ("no runinfo and an empty input dir", no_size))
    for label, skipped in skipped_samples:
        if skipped:
            print(f"  not planned, {label}: {len(skipped)} ({', '.join(skipped[:5])})", flush=True)


def write_batch_csv(path, batch):
    with open(path, "w", newline="") as f:
        csv.writer(f).writerows([("sra_id", "reserve_bytes")] + batch)


def plan(args):
    states, batches, batch, too_large, no_size = next_batch(args)
    counts = Counter(states.values())
    print(f"{args.ids}: {len(states)} samples, {len(batches)} batches in {args.staging_dir}")
    for state in STATES:
        print(f"  {state:<15} {counts[state]:>8}")
    print()
    print_batch(f"next batch {next_batch_name(args.staging_dir)}", batch, too_large, no_size, args.disk_budget)
    if args.out:
        write_batch_csv(args.out, batch)
        print(f"  written to {args.out}")


# --- results of a batch ---

def read_trace(batch_dir):
    """Return {sra_id: {process: last trace row}} from all runs of the batch."""
    steps = defaultdict(dict)
    for path in sorted(glob.glob(os.path.join(batch_dir, "pipeline_info", "*", "execution_trace.tsv"))):
        for row in read_csv(path, delimiter="\t"):
            steps[row["tag"]][row["process"]] = row
    return steps


def passes_size_gate(batch_dir, sra_id):
    published = os.path.join(batch_dir, "published")
    sylph = os.path.join(published, "branch_1", "01_taxonomy", sra_id, f"{sra_id}.sylph_profile.tsv")
    contigs = os.path.join(published, "04_shared_assembly", sra_id, "final.contigs.fa")
    return file_size(sylph) > EMPTY_SYLPH_BYTES and file_size(contigs) > 0


def sample_result(steps, gate):
    """Return (state, reason) from the last trace row of each step and the size gate."""
    failed = {process: row["exit"] for process, row in steps.items() if row["status"] not in OK_STATUS}
    permanent = {process: code for process, code in failed.items() if code not in TRANSIENT_EXITS}
    if failed:
        reason = ", ".join(f"{process} exit {code}" for process, code in sorted(failed.items()))
        return ("failed" if permanent else "retry"), reason
    missing = ENDPOINTS - set(steps)
    if missing:
        return "retry", "not run: " + ", ".join(sorted(missing))
    if not gate:
        return "failed", "size gate"
    return "done", ""


def copy_failed_logs(batch_dir, steps):
    for sra_id, rows in steps.items():
        for process, row in rows.items():
            if row["status"] in OK_STATUS:
                continue
            target = os.path.join(batch_dir, "failed_logs", process, sra_id)
            os.makedirs(target, exist_ok=True)
            for name in TASK_LOGS:
                if os.path.isfile(os.path.join(row["workdir"], name)):
                    shutil.copy2(os.path.join(row["workdir"], name), target)


def finish_batch(batch_dir, nextflow, completed):
    """Write samples.tsv, copy the logs of failed tasks and remove the work dir of the last run."""
    ids = [row["sra_id"] for row in read_csv(os.path.join(batch_dir, "batch.csv"))]
    steps = read_trace(batch_dir)
    rows = []
    for sra_id in ids:
        gate = passes_size_gate(batch_dir, sra_id)
        state, reason = sample_result(steps.get(sra_id, {}), gate)
        if state == "failed" and not completed:
            state, reason = "retry", f"run stopped; {reason}"
        rows.append((sra_id, state, "yes" if gate else "no", reason))
    copy_failed_logs(batch_dir, steps)
    temporary = os.path.join(batch_dir, "samples.tsv.tmp")
    with open(temporary, "w", newline="") as f:
        csv.writer(f, delimiter="\t").writerows([("sra_id", "state", "gate", "reason")] + rows)
    os.replace(temporary, os.path.join(batch_dir, "samples.tsv"))
    subprocess.run([nextflow, "clean", "-f", "-q"])
    counts = Counter(state for _, state, _, _ in rows)
    summary = ", ".join(f"{counts[state]} {state}" for state in ("done", "retry", "failed"))
    print(f"{os.path.basename(batch_dir)} finished: {summary}", flush=True)


def run(args):
    for batch_dir, _, states in read_batches(args.staging_dir):
        if states is None:
            print(f"{os.path.basename(batch_dir)} has no samples.tsv, so its run stopped. Finish it first.", flush=True)
            finish_batch(batch_dir, args.nextflow, completed=False)

    nextflow_run, signals = {}, []

    def pass_on(signum, frame):
        signals.append(signum)
        if "process" in nextflow_run:
            nextflow_run["process"].send_signal(signum)

    signal.signal(signal.SIGTERM, pass_on)
    signal.signal(signal.SIGINT, pass_on)
    n_run = 0
    while (args.n_batches is None or n_run < args.n_batches) and not signals:
        _, _, batch, too_large, no_size = next_batch(args)
        if not batch:
            print("no eligible sample fits the budget: stop")
            return
        name = next_batch_name(args.staging_dir)
        batch_dir = os.path.join(args.staging_dir, name)
        os.makedirs(batch_dir)
        batch_csv = os.path.join(batch_dir, "batch.csv")
        write_batch_csv(batch_csv, batch)
        print_batch(name, batch, too_large, no_size, args.disk_budget)

        nextflow_run["process"] = subprocess.Popen([
            args.nextflow, "-log", os.path.join(batch_dir, "nextflow.log"), "run", args.pipeline,
            "--input", batch_csv, "--indir", args.input_dir, "--outdir", batch_dir,
            "--limit", str(len(batch)), "-ansi-log", "false"])
        returncode = nextflow_run["process"].wait()
        if not read_trace(batch_dir):
            os.remove(batch_csv)
            raise SystemExit(f"{name}: nextflow exited with {returncode} before any task ended. "
                             f"See {batch_dir}/nextflow.log")
        finish_batch(batch_dir, args.nextflow, completed=returncode == 0 and not signals)
        n_run += 1
    if signals:
        raise SystemExit(f"stopped by signal {signals[0]}")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["plan", "run"])
    parser.add_argument("--ids", required=True, help="CSV with an sra_id column, in the order to process")
    parser.add_argument("--input-dir", required=True, help="dir with one subdir for each downloaded sample")
    parser.add_argument("--staging-dir", required=True, help="dir with the batch dirs batch_NNNN/")
    parser.add_argument("--runinfo", help="default: <input-dir>.download/runinfo.csv")
    parser.add_argument("--disk-budget", required=True, type=parse_size, help="for example 20T (T = 10^12 bytes)")
    parser.add_argument("--reserve-factor", type=float, default=5.3, help="reserve / FASTQ size from runinfo")
    parser.add_argument("--input-factor", type=float, default=36.0,
                        help="reserve / size of the input dir, when runinfo has no bases")
    parser.add_argument("--max-batches", type=int, default=2, help="batches for each sample before retry limit")
    parser.add_argument("--out", help="plan: write the next batch to this CSV")
    parser.add_argument("--n-batches", type=int, help="run: stop after this many batches (default: no limit)")
    parser.add_argument("--pipeline", default=os.path.normpath(DEFAULT_PIPELINE), help="run: the main .nf file")
    parser.add_argument("--nextflow", default="nextflow", help="run: the nextflow command")
    args = parser.parse_args()
    args.runinfo = args.runinfo or os.path.normpath(args.input_dir) + ".download/runinfo.csv"
    return args


def main():
    args = parse_args()
    if args.command == "plan":
        plan(args)
    elif args.command == "run":
        run(args)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run the pipeline in batches that fit a disk budget, and promote the results.

    batches.py plan    --staging-dir S --ids IDS.csv --input-dir IN --disk-budget 20T [--out next.csv]
    batches.py run     --staging-dir S --ids IDS.csv --input-dir IN --disk-budget 20T [--n-batches N]
    batches.py status  --staging-dir S --ids IDS.csv --input-dir IN --disk-budget 20T [--batch batch_NNNN]
    batches.py sample  --staging-dir S SRA_ID
    batches.py promote --staging-dir S --valid-dir VALID --logs-dir LOGS [--dry-run]

<staging-dir>/samples.tsv is the only record of the sample states. It has one row for
each sample that was in a batch: state, last batch, number of batches, size gate, the
sizes of the sylph profile and the contigs, reason, and one column for each pipeline
step (ok, failed <exit code>, or - when the step did not run). Only one run or promote
at a time can use a staging dir.

plan: each sample reserves --reserve-factor x its FASTQ size from runinfo
(2 x bases + 120 B x spots), or else --input-factor x the size of its dir in
--input-dir. The next batch takes the eligible samples in the order of --ids until
the next sample does not fit in --disk-budget. plan and run write the samples in
state too large or no size to <staging-dir>/not_planned.tsv.

run: creates <staging-dir>/batch_NNNN/ with batch.csv, sets its samples to running and
runs the pipeline with that dir as outdir. Then it writes the results to samples.tsv,
copies the logs of failed tasks to failed_logs/, removes from published/ the outputs
of the samples that promote does not take (their logs stay), rebuilds <staging-dir>/traces.tsv
(the traces of all batches, with the columns batch and run) and runs `nextflow clean
-f` for the session of the batch. The traces stay in the batch dirs. It repeats until
no sample is eligible or after --n-batches batches. Start it from the launch dir of the
pipeline. At its start, it finishes the batches of a stopped run. --n-batches 0 only
does that.

status: the sample states, the result of each last batch and the causes of the failed
samples. sample: the row of one sample.

promote: moves the dirs of the samples with gate yes and state done or failed from
<batch>/published/ to the same paths in --valid-dir and --logs-dir:
  <step>/<id>/logs/                              -> <logs-dir>/<step>/<id>/logs/
  <step>/<id>/ of 02_fastp_qc, 03_non_human_reads,
  or with only logs/ in it                       -> <logs-dir>/<step>/<id>/
  all other <step>/<id>/                         -> <valid-dir>/<step>/<id>/
It uses os.rename, so all dirs must be on one filesystem. If a target exists, no dir
of that sample moves. It writes each move to <staging-dir>/promote.csv and sets the
state promoted. With --dry-run, it only prints the counts. After a stop, it
moves only what is still in published/.

Sample states:
  running, done, failed, stopped, promoted   the state in samples.tsv
  not downloaded  not in samples.tsv and not in --input-dir
  too large       the reserve is larger than --disk-budget
  no size         no bases in runinfo and an empty dir in --input-dir
  eligible        new samples in --input-dir, and stopped samples

The result of a batch for each sample:
  done     no task failed and the size gate passed
  failed   a task failed, or the size gate failed. This is final. A retry is set in
           the code of a process, with more resources, for example FILTER_AND_ANNOTATE.
  stopped  the run did not end normally, and a task of the sample has the status
           ABORTED, or no task failed and an endpoint did not run. A later batch runs
           the sample again.
If the run ended normally, a step that did not run was skipped on purpose, for example
because MEGAHIT made no contigs. After a stop, such a sample runs again.
promote takes the samples with gate yes and state done or failed. run removes the
outputs of all other samples of the batch: stopped, and gate no.

Python 3.7 or later, standard library only.
"""
import argparse
import csv
import glob
import os
import re
import shutil
import subprocess
import time
from collections import Counter, defaultdict

RUN_ACCESSION = re.compile(r"[SED]RR\d+")
BATCH_NAME = re.compile(r"batch_(\d+)")
FASTQ_BYTES_PER_SPOT = 120
SIZE_UNITS = {"": 1, "K": 10**3, "M": 10**6, "G": 10**9, "T": 10**12}
STATES = ["running", "done", "failed", "promoted", "not downloaded", "too large", "no size", "eligible"]
NOT_PLANNED_STATES = ("too large", "no size")
FINAL_STATES = {"done", "failed", "promoted"}
TABLE = "samples.tsv"
TRACES = "traces.tsv"
NOT_PLANNED = "not_planned.tsv"
LOG_STEPS = {"02_fastp_qc", "03_non_human_reads"}
SAMPLE_COLUMNS = ["sra_id", "state", "batch", "batches", "gate", "sylph_bytes", "contigs_bytes", "reason"]
STEPS = [
    "CHECK_AND_CONVERT_READS", "FASTP_QC", "REMOVE_HUMAN_READS", "MEGAHIT_ASSEMBLY", "SYLPH_TAXONOMY",
    "BAKTA_ANNOTATION", "CREATE_CATALOG_AND_INDEX", "ALIGN_AND_QUANTIFY_READS", "CALCULATE_TPM_AND_ANNOTATE",
    "MAP_FOR_BINNING", "METABAT2_BINNING", "CHECKM_QA", "FILTER_AND_ANNOTATE",
]
OK_STATUS = {"COMPLETED", "CACHED"}
ENDPOINTS = {"SYLPH_TAXONOMY", "CALCULATE_TPM_AND_ANNOTATE", "FILTER_AND_ANNOTATE"}
# A sylph profile with only its header has 218 bytes.
EMPTY_SYLPH_BYTES = 218
TASK_LOGS = [".command.sh", ".command.err", ".command.out", ".command.log", ".exitcode"]
DEFAULT_PIPELINE = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "pipeline", "mtg_pipeline.nf")


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


# --- the status table and sample states ---

def read_table(staging_dir):
    """Return {sra_id: row} of <staging-dir>/samples.tsv."""
    path = os.path.join(staging_dir, TABLE)
    return {row["sra_id"]: row for row in read_csv(path, delimiter="\t")} if os.path.exists(path) else {}


def write_table(staging_dir, table):
    other_steps = sorted({key for row in table.values() for key in row} - set(SAMPLE_COLUMNS) - set(STEPS))
    path = os.path.join(staging_dir, TABLE)
    with open(path + ".tmp", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=SAMPLE_COLUMNS + STEPS + other_steps, delimiter="\t", restval="-")
        writer.writeheader()
        writer.writerows(table.values())
    os.replace(path + ".tmp", path)


def next_batch_name(staging_dir):
    numbers = [0]
    for name in subdirs(staging_dir):
        match = BATCH_NAME.fullmatch(name)
        if match:
            numbers.append(int(match.group(1)))
    return f"batch_{max(numbers) + 1:04d}"


def sample_states(ids, table, input_dir):
    """Return {sra_id: one of STATES}."""
    downloaded = subdirs(input_dir)
    result = {}
    for sra_id in ids:
        row = table.get(sra_id)
        if row and row["state"] in FINAL_STATES | {"running"}:
            result[sra_id] = row["state"]
        elif sra_id not in downloaded:
            result[sra_id] = "not downloaded"
        else:
            result[sra_id] = "eligible"
    return result


def reserves(args, states):
    """Return {sra_id: reserved bytes} for the eligible samples."""
    fastq_sizes = read_fastq_sizes(args.runinfo or os.path.normpath(args.input_dir) + ".download/runinfo.csv")

    def reserve_of(sra_id):
        if sra_id in fastq_sizes:
            return int(args.reserve_factor * fastq_sizes[sra_id])
        return int(args.input_factor * dir_bytes(os.path.join(args.input_dir, sra_id)))

    return {sra_id: reserve_of(sra_id) for sra_id, state in states.items() if state == "eligible"}


def load_states(args):
    """Return (status table, {sra_id: state} for the samples of --ids, {sra_id: reserved bytes})."""
    table = read_table(args.staging_dir)
    states = sample_states(read_ids(args.ids), table, args.input_dir)
    reserve = reserves(args, states)
    for sra_id, reserved in reserve.items():
        if reserved == 0:
            states[sra_id] = "no size"
        elif reserved > args.disk_budget:
            states[sra_id] = "too large"
    return table, states, reserve


def choose_batch(states, reserve, budget):
    """Return [(sra_id, reserved bytes)] for the next batch."""
    batch, total = [], 0
    for sra_id, state in states.items():
        if state != "eligible":
            continue
        if total + reserve[sra_id] > budget:
            break
        batch.append((sra_id, reserve[sra_id]))
        total += reserve[sra_id]
    return batch


def write_not_planned(staging_dir, states, reserve):
    """Write the samples in state too large or no size to <staging-dir>/not_planned.tsv."""
    rows = [(sra_id, state, reserve[sra_id]) for sra_id, state in states.items() if state in NOT_PLANNED_STATES]
    os.makedirs(staging_dir, exist_ok=True)
    path = os.path.join(staging_dir, NOT_PLANNED)
    with open(path, "w", newline="") as f:
        csv.writer(f, delimiter="\t").writerows([("sra_id", "state", "reserve_bytes")] + rows)
    if rows:
        counts = Counter(state for _, state, _ in rows)
        print(f"  not planned: {counts['too large']} too large, {counts['no size']} no size, see {path}", flush=True)


def print_batch(name, batch, budget):
    reserved = sum(size for _, size in batch)
    print(f"{name}: {len(batch)} samples, reserve {reserved / 1e9:,.1f} GB of {budget / 1e9:,.1f} GB", flush=True)


def write_batch_csv(path, batch):
    with open(path, "w", newline="") as f:
        csv.writer(f).writerows([("sra_id", "reserve_bytes")] + batch)


def print_states(args, states, table):
    counts = Counter(states.values())
    print(f"{args.ids}: {len(states)} samples, {len(table)} of them in {os.path.join(args.staging_dir, TABLE)}")
    for state in STATES:
        print(f"  {state:<15} {counts[state]:>8}")
    print()


def plan(args):
    table, states, reserve = load_states(args)
    print_states(args, states, table)
    print_batch(f"next batch {next_batch_name(args.staging_dir)}", choose_batch(states, reserve, args.disk_budget),
                args.disk_budget)
    write_not_planned(args.staging_dir, states, reserve)
    if args.out:
        write_batch_csv(args.out, choose_batch(states, reserve, args.disk_budget))
        print(f"  written to {args.out}")


def show_sample(args):
    row = read_table(args.staging_dir).get(args.sra_id, {"sra_id": args.sra_id, "state": "not in any batch"})
    for column, value in row.items():
        print(f"{column:<28} {value}")


def status(args):
    table, states, _ = load_states(args)
    print_states(args, states, table)
    rows = [row for row in table.values() if not args.batch or row["batch"] == args.batch]

    shown = ("running", "done", "stopped", "failed", "promoted")
    print(f"{'last batch':<12} {'samples':>8}" + "".join(f" {state:>8}" for state in shown) + f" {'gate yes':>9}")
    for name in sorted({row["batch"] for row in rows}):
        batch_rows = [row for row in rows if row["batch"] == name]
        counts = Counter(row["state"] for row in batch_rows)
        gate = sum(row["gate"] == "yes" for row in batch_rows)
        print(f"{name:<12} {len(batch_rows):>8}" + "".join(f" {counts[state]:>8}" for state in shown) + f" {gate:>9}")

    causes = Counter(cause for row in rows if row["state"] == "failed" for cause in row["reason"].split(", "))
    print("\nfailed samples by cause (one sample can have more than one cause):")
    for cause, count in causes.most_common():
        print(f"  {count:>8}  {cause}")


# --- results of a batch ---

def read_trace(batch_dir):
    """Return {sra_id: {process: last trace row}} from all runs of the batch."""
    steps = defaultdict(dict)
    for path in sorted(glob.glob(os.path.join(batch_dir, "pipeline_info", "*", "execution_trace.tsv"))):
        for row in read_csv(path, delimiter="\t"):
            steps[row["tag"]][row["process"]] = row
    return steps


def gate_sizes(published, sra_id):
    """Return (bytes of the sylph profile, bytes of the contigs) of one sample in published/."""
    sylph = os.path.join(published, "branch_1", "01_taxonomy", sra_id, f"{sra_id}.sylph_profile.tsv")
    contigs = os.path.join(published, "04_shared_assembly", sra_id, "final.contigs.fa")
    return file_size(sylph), file_size(contigs)


def sample_result(steps, gate, completed):
    """Return (state, reason) from the last trace row of each step, the size gate and how the run ended."""
    failed = {process: row["exit"] for process, row in steps.items() if row["status"] not in OK_STATUS}
    aborted = any(row["status"] == "ABORTED" for row in steps.values())
    if not completed and (aborted or (not failed and not ENDPOINTS <= set(steps))):
        return "stopped", "run stopped"
    if failed:
        return "failed", ", ".join(f"{process} exit {code}" for process, code in sorted(failed.items()))
    if not gate:
        return "failed", "size gate"
    return "done", ""


def step_status(trace_row):
    return "ok" if trace_row["status"] in OK_STATUS else f"failed {trace_row['exit']}"


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


def write_traces(staging_dir):
    """Rebuild <staging-dir>/traces.tsv from the traces of all batches."""
    paths = sorted(glob.glob(os.path.join(staging_dir, "batch_*", "pipeline_info", "*", "execution_trace.tsv")))
    columns = ["batch", "run"]
    for path in paths:
        with open(path, newline="") as f:
            columns += [column for column in next(csv.reader(f, delimiter="\t"), []) if column not in columns]
    target = os.path.join(staging_dir, TRACES)
    with open(target + ".tmp", "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for path in paths:
            run_dir = os.path.dirname(path)
            batch = os.path.basename(os.path.dirname(os.path.dirname(run_dir)))
            with open(path, newline="") as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    writer.writerow({"batch": batch, "run": os.path.basename(run_dir), **row})
    os.replace(target + ".tmp", target)


def start_batch(staging_dir, table, name, ids):
    """Write a copy of the table with the samples of the batch in state running."""
    started = dict(table)
    for sra_id in ids:
        batches = int(table.get(sra_id, {}).get("batches", 0)) + 1
        started[sra_id] = {"sra_id": sra_id, "state": "running", "batch": name, "batches": batches}
    write_table(staging_dir, started)


def remove_outputs(published, steps, sra_id):
    """Remove the outputs of one sample from published/. Keep its logs."""
    for step in steps:
        sample_dir = os.path.join(published, step, sra_id)
        if step in LOG_STEPS or not os.path.isdir(sample_dir):
            continue
        for name in os.listdir(sample_dir):
            path = os.path.join(sample_dir, name)
            if name == "logs":
                continue
            if os.path.isdir(path) and not os.path.islink(path):
                shutil.rmtree(path)
            else:
                os.remove(path)


def session_id(batch_dir):
    """Return the Nextflow session id of the batch run from its nextflow.log, or None."""
    with open(os.path.join(batch_dir, "nextflow.log")) as f:
        for line in f:
            if "Session UUID: " in line:
                return line.split("Session UUID: ")[1].strip()
    return None


def finish_batch(staging_dir, name, nextflow, completed):
    """Copy the logs of failed tasks, remove the outputs that promote does not take, remove the work dir and
    write the results of the batch to samples.tsv. samples.tsv comes last: after a kill in this function, the
    next run does all steps again."""
    batch_dir = os.path.join(staging_dir, name)
    published = os.path.join(batch_dir, "published")
    ids = [row["sra_id"] for row in read_csv(os.path.join(batch_dir, "batch.csv"))]
    steps = read_trace(batch_dir)
    table = read_table(staging_dir)
    counts = Counter()
    not_promoted = []
    for sra_id in ids:
        sylph, contigs = gate_sizes(published, sra_id)
        gate = sylph > EMPTY_SYLPH_BYTES and contigs > 0
        state, reason = sample_result(steps.get(sra_id, {}), gate, completed)
        table[sra_id] = {
            "sra_id": sra_id, "state": state, "batch": name, "batches": table.get(sra_id, {}).get("batches", 1),
            "gate": "yes" if gate else "no", "sylph_bytes": sylph, "contigs_bytes": contigs, "reason": reason,
            **{process: step_status(row) for process, row in steps.get(sra_id, {}).items()},
        }
        counts[state] += 1
        if state == "stopped" or not gate:
            not_promoted.append(sra_id)
    copy_failed_logs(batch_dir, steps)
    published_steps = step_dirs(published)
    for sra_id in not_promoted:
        remove_outputs(published, published_steps, sra_id)
    write_traces(staging_dir)
    session = session_id(batch_dir)
    if session:
        subprocess.run([nextflow, "clean", "-f", "-q", session], check=True)
    write_table(staging_dir, table)
    summary = ", ".join(f"{counts[state]} {state}" for state in ("done", "stopped", "failed"))
    print(f"{name} finished: {summary}; outputs removed for {len(not_promoted)} samples", flush=True)


def run(args):
    for name in sorted({row["batch"] for row in read_table(args.staging_dir).values() if row["state"] == "running"}):
        print(f"{name} has samples in state running, so its run stopped. Finish it first.", flush=True)
        finish_batch(args.staging_dir, name, args.nextflow, completed=False)

    n_run = 0
    while args.n_batches is None or n_run < args.n_batches:
        table, states, reserve = load_states(args)
        batch = choose_batch(states, reserve, args.disk_budget)
        if not batch:
            print("no eligible sample: stop", flush=True)
            write_not_planned(args.staging_dir, states, reserve)
            return
        name = next_batch_name(args.staging_dir)
        batch_dir = os.path.join(args.staging_dir, name)
        os.makedirs(batch_dir)
        batch_csv = os.path.join(batch_dir, "batch.csv")
        write_batch_csv(batch_csv, batch)
        start_batch(args.staging_dir, table, name, [sra_id for sra_id, _ in batch])
        print_batch(name, batch, args.disk_budget)
        write_not_planned(args.staging_dir, states, reserve)

        # Not subprocess.run: on Ctrl-C it kills nextflow before nextflow can cancel its tasks.
        returncode = subprocess.Popen([
            args.nextflow, "-log", os.path.join(batch_dir, "nextflow.log"), "run", args.pipeline,
            "--input", batch_csv, "--indir", args.input_dir, "--outdir", batch_dir,
            "--limit", str(len(batch)), "-ansi-log", "false"]).wait()
        if not read_trace(batch_dir):
            write_table(args.staging_dir, table)
            raise SystemExit(f"{name}: nextflow exited with {returncode} before any task ended. "
                             f"See {batch_dir}/nextflow.log")
        finish_batch(args.staging_dir, name, args.nextflow, completed=returncode == 0)
        n_run += 1


# --- promote ---

def step_dirs(published):
    """Return the step dirs of published/ relative to it, for example branch_3/02_binning."""
    steps = []
    for name in sorted(subdirs(published)):
        if name.startswith("branch_"):
            steps += [os.path.join(name, step) for step in sorted(subdirs(os.path.join(published, name)))]
        else:
            steps.append(name)
    return steps


def sample_moves(published, steps, sra_id, valid_dir, logs_dir):
    """Return [(source, target)] for the dirs of one sample that are still in published/."""
    moves = []
    for step in steps:
        source = os.path.join(published, step, sra_id)
        if not os.path.isdir(source):
            continue
        if step in LOG_STEPS or set(os.listdir(source)) <= {"logs"}:
            moves.append((source, os.path.join(logs_dir, step, sra_id)))
            continue
        logs = os.path.join(source, "logs")
        if os.path.isdir(logs):
            moves.append((logs, os.path.join(logs_dir, step, sra_id, "logs")))
        moves.append((source, os.path.join(valid_dir, step, sra_id)))
    return moves


def plan_promote(args, table):
    """Return ([(row, moves)] ready to promote, conflicts, samples with nothing to move)."""
    ready, conflicts, missing, steps_of = [], [], [], {}
    for row in table.values():
        if row["gate"] != "yes" or row["state"] not in ("done", "failed"):
            continue
        published = os.path.join(args.staging_dir, row["batch"], "published")
        if published not in steps_of:
            steps_of[published] = step_dirs(published)
        moves = sample_moves(published, steps_of[published], row["sra_id"], args.valid_dir, args.logs_dir)
        taken = [target for _, target in moves if os.path.exists(target)]
        if taken:
            conflicts.append(f"{row['sra_id']}: {taken[0]} exists")
        elif moves or os.path.isdir(os.path.join(args.valid_dir, "04_shared_assembly", row["sra_id"])):
            ready.append((row, moves))
        else:
            missing.append(row["sra_id"])
    return ready, conflicts, missing


def promote(args):
    table = read_table(args.staging_dir)
    if any(row["state"] == "running" for row in table.values()):
        raise SystemExit("samples.tsv has samples in state running. Wait until run ends. "
                         "After a stop, run --n-batches 0 finishes the stopped batch.")
    ready, conflicts, missing = plan_promote(args, table)
    targets = [target for _, moves in ready for _, target in moves]
    to_valid = sum(target.startswith(os.path.join(args.valid_dir, "")) for target in targets)
    states = Counter(row["state"] for row, _ in ready)
    print(f"ready      {len(ready):>8}  (done {states['done']}, failed {states['failed']}), "
          f"{to_valid} dirs to {args.valid_dir}, {len(targets) - to_valid} to {args.logs_dir}")
    print(f"conflict   {len(conflicts):>8}  a target exists, so nothing of the sample moves")
    print(f"missing    {len(missing):>8}  nothing in published/ and no 04_shared_assembly in --valid-dir")
    for line in conflicts[:5] + missing[:5]:
        print(f"  {line}")
    if args.dry_run:
        print("dry run: nothing moved.")
        return

    log_path = os.path.join(args.staging_dir, "promote.csv")
    new_log = not os.path.exists(log_path)
    with open(log_path, "a", newline="", buffering=1) as f:
        log = csv.writer(f)
        if new_log:
            log.writerow(["time_utc", "sra_id", "batch", "source", "target"])
        try:
            for n, (row, moves) in enumerate(ready, 1):
                for source, target in moves:
                    os.makedirs(os.path.dirname(target), exist_ok=True)
                    os.rename(source, target)
                    log.writerow([time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), row["sra_id"], row["batch"],
                                  source, target])
                row["state"] = "promoted"
                if n % 1000 == 0:
                    write_table(args.staging_dir, table)
        finally:
            write_table(args.staging_dir, table)
    print(f"promoted {len(ready)} samples")


def parse_args():
    staging = argparse.ArgumentParser(add_help=False)
    staging.add_argument("--staging-dir", required=True, help="dir with samples.tsv and the batch dirs batch_NNNN/")
    samples = argparse.ArgumentParser(add_help=False)
    samples.add_argument("--ids", required=True, help="CSV with an sra_id column, in the order to process")
    samples.add_argument("--input-dir", required=True, help="dir with one subdir for each downloaded sample")
    budget = argparse.ArgumentParser(add_help=False)
    budget.add_argument("--disk-budget", type=parse_size, required=True, help="for example 20T (T = 10^12 bytes)")
    budget.add_argument("--runinfo", help="default: <input-dir>.download/runinfo.csv")
    budget.add_argument("--reserve-factor", type=float, default=4.5, help="reserve / FASTQ size from runinfo")
    budget.add_argument("--input-factor", type=float, default=36.0,
                        help="reserve / size of the input dir, when runinfo has no bases")

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("plan", parents=[staging, samples, budget], help="print the next batch")
    command.add_argument("--out", help="write the next batch to this CSV")
    command.set_defaults(func=plan)
    command = commands.add_parser("run", parents=[staging, samples, budget], help="run batches")
    command.add_argument("--n-batches", type=int, help="stop after this many batches (default: no limit)")
    command.add_argument("--pipeline", default=os.path.normpath(DEFAULT_PIPELINE), help="the main .nf file")
    command.add_argument("--nextflow", default="nextflow", help="the nextflow command")
    command.set_defaults(func=run)
    command = commands.add_parser("status", parents=[staging, samples, budget], help="print the sample states")
    command.add_argument("--batch", help="only the samples whose last batch is this one")
    command.set_defaults(func=status)
    command = commands.add_parser("sample", parents=[staging], help="print the row of one sample")
    command.add_argument("sra_id")
    command.set_defaults(func=show_sample)
    command = commands.add_parser("promote", parents=[staging], help="move the accepted samples")
    command.add_argument("--valid-dir", required=True, help="the target of the results, for example gutmap/valid")
    command.add_argument("--logs-dir", required=True, help="the target of the logs, for example gutmap/valid-logs")
    command.add_argument("--dry-run", action="store_true", help="only print the counts, move nothing")
    command.set_defaults(func=promote)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    args.func(args)

#!/usr/bin/env python3
"""Download SRA/ENA runs with kingfisher into the input dir that the pipeline reads.

    download.py runinfo  --ids samples.csv --input-dir INPUT_DIR
    download.py plan     --ids samples.csv --input-dir INPUT_DIR
    download.py run      --ids samples.csv --input-dir INPUT_DIR [--methods aws-http:4,ena-ftp:4]
    download.py validate --input-dir INPUT_DIR [--threads 4]
    download.py status   --input-dir INPUT_DIR [--job JOBID]

The state dir (default: <input-dir>.download) holds the attempt records, the NCBI
runinfo of the samples and the downloads that are not in the input dir yet. runinfo
adds the rows of new samples to <state-dir>/runinfo.csv. run puts each download that passes the
file-set check in <state-dir>/downloaded/. validate reads all data of these downloads
and moves each one that passes to the input dir.
"""
import argparse
import csv
import os
import re
import shutil
import signal
import socket
import subprocess
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed

import runinfo
import validate

RUN_ACCESSION = re.compile(r"[SED]RR\d+")
COUNTED_RESULTS = {"failed", "invalid", "timeout"}
ATTEMPT_RESULTS = ["ok", "failed", "invalid", "timeout", "interrupted"]
ATTEMPT_COLUMNS = ["sra_id", "attempt", "job_id", "method", "started_utc", "seconds",
                   "result", "bytes", "files", "reason"]
VALIDATION_COLUMNS = ["sra_id", "job_id", "started_utc", "seconds", "result", "bytes", "reason"]
CAUSE_PATTERNS = [
    r"errorCode=\d+ (?!URI=)[^\n]+",                          # aria2c
    r"ERROR: (?!No more specified download methods)[^\n]+",   # kingfisher
]
PROGRESS_SECONDS = 60
RUNINFO_BATCH = 200


def utc_now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def job_id():
    return os.environ.get("SLURM_JOB_ID") or f"local-{socket.gethostname()}-{os.getpid()}"


# --- input and records ---

def expected_bytes(value):
    try:
        size = int(float(value))
    except (TypeError, ValueError):
        return None
    return size if size > 0 else None


def read_samples(path):
    """Return {sra_id: expected bytes or None} and the counts of the rows that are not used."""
    samples, empty_rows, duplicates, not_runs = {}, 0, 0, 0
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        if "sra_id" not in (reader.fieldnames or []):
            raise SystemExit(f"{path}: no sra_id column")
        for row in reader:
            sra_id = (row["sra_id"] or "").strip()
            if not sra_id:
                empty_rows += 1
            elif not RUN_ACCESSION.fullmatch(sra_id):
                not_runs += 1
            elif sra_id in samples:
                duplicates += 1
            else:
                samples[sra_id] = expected_bytes(row.get("Bytes"))
    return samples, empty_rows, duplicates, not_runs


def read_records(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def append_record(path, columns, record):
    is_new = not os.path.exists(path)
    with open(path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        if is_new:
            writer.writeheader()
        writer.writerow(record)


def read_runinfo(state_dir):
    """Return {sra_id: runinfo row} from <state-dir>/runinfo.csv."""
    return {row["Run"]: row for row in read_records(os.path.join(state_dir, "runinfo.csv"))}


def subdirs(path):
    return set(os.listdir(path)) if os.path.isdir(path) else set()


def dir_bytes(path):
    return sum(os.path.getsize(os.path.join(path, name)) for name in os.listdir(path))


def read_history(state_dir):
    """Return the rows of attempts.csv and the counted failures of each sample.

    The counted failures are failed download attempts plus failed content checks.
    """
    attempt_rows = read_records(os.path.join(state_dir, "attempts.csv"))
    validation_rows = read_records(os.path.join(state_dir, "validations.csv"))
    failures = Counter(row["sra_id"] for row in attempt_rows if row["result"] in COUNTED_RESULTS)
    failures.update(row["sra_id"] for row in validation_rows if row["result"] != "ok")
    return attempt_rows, failures


def classify(samples, input_dir, state_dir, failures, max_attempts):
    # downloaded/ first: validate can move a dir to the input dir between the two listings.
    downloaded = subdirs(os.path.join(state_dir, "downloaded"))
    present = subdirs(input_dir)
    classes = {"present": [], "downloaded": [], "exhausted": [], "to download": []}
    for sra_id in samples:
        if sra_id in present:
            classes["present"].append(sra_id)
        elif sra_id in downloaded:
            classes["downloaded"].append(sra_id)
        elif failures[sra_id] >= max_attempts:
            classes["exhausted"].append(sra_id)
        else:
            classes["to download"].append(sra_id)
    return classes


def print_class_row(name, ids, samples, note=""):
    gb = sum(samples[sra_id] or 0 for sra_id in ids) / 1e9
    print(f"{name:<12} {len(ids):>9} {gb:>10.1f}  {note}".rstrip())


def plan(args, failures):
    """Print the plan and return the sample ids of each class."""
    samples, empty_rows, duplicates, not_runs = read_samples(args.ids)
    from_csv = sum(size is not None for size in samples.values())
    runinfo_rows = read_runinfo(args.state_dir)
    for sra_id, size in samples.items():
        megabytes = expected_bytes(runinfo_rows.get(sra_id, {}).get("size_MB"))
        if size is None and megabytes:
            samples[sra_id] = megabytes * 10**6
    from_runinfo = sum(size is not None for size in samples.values()) - from_csv
    classes = classify(samples, args.input_dir, args.state_dir, failures, args.max_attempts)
    n_rows = len(samples) + empty_rows + duplicates + not_runs
    notes = {"downloaded": "wait for the content check",
             "exhausted": f"{args.max_attempts} or more counted attempts"}

    print(f"{args.ids}: {len(samples)} samples from {n_rows} rows "
          f"(not used: {empty_rows} empty, {duplicates} duplicate, {not_runs} not a run accession)")
    print(f"sizes: {from_csv} from the Bytes column, {from_runinfo} from runinfo.csv, "
          f"{len(samples) - from_csv - from_runinfo} unknown")
    print(f"\n{'class':<12} {'samples':>9} {'GB':>10}")
    for name, ids in classes.items():
        print_class_row(name, ids, samples, notes.get(name, ""))
    print_class_row("total", list(samples), samples)
    return classes


# --- runinfo ---

def fetch_batch(batch, credentials, pacer):
    """Return (runinfo rows, {sra_id: error}). If the request for the batch keeps failing, ask for each sample alone.

    NCBI fails on some samples every time (SRR769534: HTTP 500 alone, a cut response in a batch).
    """
    try:
        return runinfo.parse(runinfo.fetch(batch, credentials, pacer)), {}
    except RuntimeError as error:
        print(f"runinfo: request for {len(batch)} samples failed ({error}), asking for each sample alone", flush=True)
    rows, errors = [], {}
    for sra_id in batch:
        try:
            rows += runinfo.parse(runinfo.fetch([sra_id], credentials, pacer, retries=1))
        except RuntimeError as error:
            errors[sra_id] = str(error)
    if len(errors) == len(batch):
        raise RuntimeError(f"all {len(batch)} samples failed alone, for example {next(iter(errors.values()))}")
    return rows, errors


def fetch_runinfo(args):
    """Add the NCBI runinfo rows of the samples that are not in <state-dir>/runinfo.csv yet."""
    samples = read_samples(args.ids)[0]
    path = os.path.join(args.state_dir, "runinfo.csv")
    known = read_runinfo(args.state_dir)
    to_fetch = [sra_id for sra_id in samples if sra_id not in known]
    credentials = {name: os.environ.get(variable) for name, variable in
                   (("api_key", "NCBI_API_KEY"), ("email", "NCBI_EMAIL"), ("tool", "NCBI_TOOL"))}
    pacer = runinfo.Pacer(bool(credentials["api_key"]))
    print(f"runinfo: {len(samples)} samples, {len(samples) - len(to_fetch)} already in {path}, "
          f"{len(to_fetch)} to fetch ({'with' if credentials['api_key'] else 'without'} NCBI_API_KEY)", flush=True)
    os.makedirs(args.state_dir, exist_ok=True)

    returned, errors = set(), {}
    for start in range(0, len(to_fetch), RUNINFO_BATCH):
        batch = to_fetch[start:start + RUNINFO_BATCH]
        try:
            rows, batch_errors = fetch_batch(batch, credentials, pacer)
        except (RuntimeError, ValueError) as error:
            raise SystemExit(f"runinfo stopped after {start} samples: {error}. Run it again to continue.")
        new_rows = [row for row in rows if row[0] in batch and row[0] not in returned]
        is_new = not os.path.exists(path)
        with open(path, "a", newline="") as f:
            writer = csv.writer(f)
            if is_new:
                writer.writerow(runinfo.COLUMNS)
            writer.writerows(new_rows)
        returned.update(row[0] for row in new_rows)
        errors.update(batch_errors)
        if (start // RUNINFO_BATCH + 1) % 10 == 0:
            print(f"runinfo: {start + len(batch)}/{len(to_fetch)} samples asked, {len(returned)} rows", flush=True)

    missing = [[sra_id, errors.get(sra_id, "not returned")] for sra_id in to_fetch if sra_id not in returned]
    with open(os.path.join(args.state_dir, "runinfo_missing.csv"), "w", newline="") as f:
        csv.writer(f).writerows([["sra_id", "reason"]] + missing)
    print(f"runinfo: {len(returned)} rows added, {len(missing)} samples without runinfo "
          f"({len(errors)} of them failed) in runinfo_missing.csv")


# --- download ---

def parse_methods(text):
    """'aws-http:4,ena-ftp:4' -> [('aws-http', 4), ('ena-ftp', 4)]"""
    methods = []
    for item in text.split(","):
        name, _, weight = item.strip().partition(":")
        if not name or not weight.isdigit() or int(weight) < 1:
            raise SystemExit(f"bad --methods item {item!r}, expected name:threads")
        methods.append((name, int(weight)))
    return methods


def failure_reason(log_path):
    """Return the first known error cause in the kingfisher log, else its last line."""
    with open(log_path, errors="replace") as f:
        text = f.read()
    for pattern in CAUSE_PATTERNS:
        match = re.search(pattern, text)
        if match:
            return match.group(0).strip()[:300]
    lines = text.strip().splitlines()
    return lines[-1].strip()[:300] if lines else ""


def kill_group(process):
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass


class DownloadRun:
    """The shared state of one run: the sample queue, the kingfisher processes and the stop flag.

    In one run, a sample gets at most one attempt with each method. A failed sample goes
    back to the queue while it has a method left and fewer than max_attempts counted failures.
    """

    def __init__(self, sra_ids, args, attempt_rows, failures):
        self.args = args
        self.attempts_csv = os.path.join(args.state_dir, "attempts.csv")
        self.methods = {name for name, _ in args.methods}
        self.waiting = list(sra_ids)
        self.tried = {sra_id: set() for sra_id in sra_ids}
        self.attempts = Counter(row["sra_id"] for row in attempt_rows)
        self.failures = failures
        self.running = 0
        self.processes = set()
        self.stopping = False
        self.records = []
        self.start = time.time()
        # RLock: the SIGTERM handler runs stop() in the main thread, which can hold the lock in progress().
        self.lock = threading.RLock()

    def take(self, method):
        """Return (sra_id, attempt), or None when this method has no more work in this run."""
        while True:
            with self.lock:
                if self.stopping:
                    return None
                for sra_id in self.waiting:
                    if method not in self.tried[sra_id]:
                        self.waiting.remove(sra_id)
                        self.tried[sra_id].add(method)
                        self.attempts[sra_id] += 1
                        self.running += 1
                        return sra_id, self.attempts[sra_id]
                if self.running == 0:
                    return None
            # A running attempt can fail and put a sample back that this method has not tried.
            time.sleep(5)

    def run_kingfisher(self, cmd, cwd, log):
        """Return the exit code, or None after the attempt timeout.

        kingfisher runs in its own process group, so that the timeout and stop() also end
        its aria2c or curl children.
        """
        with self.lock:
            if self.stopping:
                return -signal.SIGTERM  # as if kingfisher got the SIGTERM
            process = subprocess.Popen(cmd, cwd=cwd, stdout=log, stderr=subprocess.STDOUT,
                                       start_new_session=True)
            self.processes.add(process)
        try:
            return process.wait(timeout=self.args.attempt_timeout)
        except subprocess.TimeoutExpired:
            kill_group(process)
            process.wait()
            return None
        finally:
            with self.lock:
                self.processes.discard(process)

    def stop(self):
        """Give out no more samples, and end the kingfisher process groups.

        The signal handler calls it: scancel, the time limit, Ctrl-C or a closed terminal.
        """
        with self.lock:
            self.stopping = True
            for process in self.processes:
                kill_group(process)

    def finish(self, record):
        sra_id = record["sra_id"]
        with self.lock:
            self.running -= 1
            self.records.append(record)
            if record["result"] in COUNTED_RESULTS:
                self.failures[sra_id] += 1
                untried = self.methods - self.tried[sra_id]
                if untried and self.failures[sra_id] < self.args.max_attempts:
                    self.waiting.append(sra_id)
            append_record(self.attempts_csv, ATTEMPT_COLUMNS, record)
            print(f"{sra_id}  {record['method']}  {record['result']}  {record['seconds']} s  "
                  f"{record['bytes'] / 1e6:.1f} MB  {record['reason']}", flush=True)

    def progress(self):
        with self.lock:
            hours = (time.time() - self.start) / 3600
            done = len(self.tried) - len(self.waiting) - self.running
            ok = [record for record in self.records if record["result"] == "ok"]
            gb = sum(record["bytes"] for record in ok) / 1e9
            remaining = len(self.tried) - done
            rate = done / hours
            eta = f"{remaining / rate:.1f} h" if done else "-"
            attempts = Counter(record["method"] for record in self.records)
            oks = Counter(record["method"] for record in ok)
            methods = ", ".join(f"{method} {oks[method]}/{attempts[method]} ok"
                                for method in sorted(self.methods))
            return (f"progress {hours:.2f} h: {done}/{len(self.tried)} samples done, {len(ok)} ok, "
                    f"{gb:.1f} GB, {gb / hours:.1f} GB/h, ETA {eta} ({methods})")


def run_attempt(sra_id, method, attempt, run_state):
    """Download one sample into staging/, check it, and move it to downloaded/. Return the record."""
    state_dir = run_state.args.state_dir
    staging = os.path.join(state_dir, "staging", sra_id)
    log_path = os.path.join(state_dir, "logs", f"{sra_id}.{attempt}.log")
    # No --check-md5sums: kingfisher 0.5.0 never compares the md5 (wwood/kingfisher-download PR #59).
    cmd = ["kingfisher", "get", "-r", sra_id, "-m", method, "-f", "fastq.gz", "sra",
           "--output-directory", staging]
    if method == "aws-http":
        cmd.append("--guess-aws-location")

    started_utc, start = utc_now(), time.time()
    os.makedirs(staging)
    with open(log_path, "w") as log:
        exit_code = run_state.run_kingfisher(cmd, staging, log)
    names = sorted(os.listdir(staging))
    size = dir_bytes(staging)

    if exit_code is None:
        result, reason = "timeout", f"no exit after {run_state.args.attempt_timeout} s"
    elif exit_code != 0 and run_state.stopping:
        result, reason = "interrupted", "the run was stopped"
    elif exit_code != 0:
        result, reason = "failed", f"kingfisher exit {exit_code}: {failure_reason(log_path)}"
    else:
        problems = validate.find_problems(staging)
        if problems:
            result, reason = "invalid", "; ".join(problems)
        else:
            result, reason = "ok", ""
    if result == "ok":
        os.rename(staging, os.path.join(state_dir, "downloaded", sra_id))
        os.remove(log_path)
    shutil.rmtree(staging, ignore_errors=True)

    return {"sra_id": sra_id, "attempt": attempt, "job_id": job_id(), "method": method,
            "started_utc": started_utc, "seconds": round(time.time() - start),
            "result": result, "bytes": size, "files": ";".join(names), "reason": reason}


def worker(method, run_state):
    while True:
        taken = run_state.take(method)
        if taken is None:
            return
        sra_id, attempt = taken
        try:
            record = run_attempt(sra_id, method, attempt, run_state)
        except OSError as exc:
            record = {"sra_id": sra_id, "attempt": attempt, "job_id": job_id(), "method": method,
                      "started_utc": utc_now(), "seconds": 0, "result": "failed", "bytes": 0,
                      "reason": f"{type(exc).__name__}: {exc}"}
        run_state.finish(record)


def check_same_filesystem(input_dir, state_dir):
    if os.stat(input_dir).st_dev != os.stat(state_dir).st_dev:
        raise SystemExit(f"{input_dir} and {state_dir} must be on one filesystem, so that a rename is atomic")


def run(args):
    """Download the samples of the plan, with one thread for each method slot.

    Only one run at a time can use a state dir. download.sbatch makes sure of this
    with --dependency=singleton.
    """
    os.makedirs(args.state_dir, exist_ok=True)
    check_same_filesystem(args.input_dir, args.state_dir)
    shutil.rmtree(os.path.join(args.state_dir, "staging"), ignore_errors=True)
    for name in ("staging", "downloaded", "logs"):
        os.makedirs(os.path.join(args.state_dir, name), exist_ok=True)

    attempt_rows, failures = read_history(args.state_dir)
    to_download = plan(args, failures)["to download"]
    run_state = DownloadRun(to_download, args, attempt_rows, failures)
    for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(signum, lambda signum, frame: run_state.stop())

    threads = [threading.Thread(target=worker, args=(name, run_state))
               for name, count in args.methods for _ in range(count)]
    for thread in threads:
        thread.start()
    next_progress = time.time() + PROGRESS_SECONDS
    while any(thread.is_alive() for thread in threads):
        time.sleep(1)
        if time.time() >= next_progress:
            print(run_state.progress(), flush=True)
            next_progress += PROGRESS_SECONDS
    print(run_state.progress(), flush=True)


# --- content check ---

def check_download(sra_id, args):
    """Read all data of one download, then move it to the input dir or remove it. Return the record."""
    sample_dir = os.path.join(args.state_dir, "downloaded", sra_id)
    started_utc, start = utc_now(), time.time()
    size = dir_bytes(sample_dir)
    problems = validate.find_problems(sample_dir, content=True, timeout=args.file_timeout)
    if problems:
        shutil.rmtree(sample_dir)
        result, reason = "invalid", "; ".join(problems)
    else:
        os.rename(sample_dir, os.path.join(args.input_dir, sra_id))
        result, reason = "ok", ""
    return {"sra_id": sra_id, "job_id": job_id(), "started_utc": started_utc,
            "seconds": round(time.time() - start), "result": result, "bytes": size, "reason": reason}


def check_downloads(args):
    """Do the content check of each dir in downloaded/, with --threads checks at the same time.

    Only one validate at a time can use a state dir. validate.sbatch makes sure of this
    with --dependency=singleton. A stop at any time leaves each dir in downloaded/ or in
    the input dir, because the check only reads and the move to the input dir is one rename.
    """
    for tool in ("vdb-validate", "gzip"):
        if not shutil.which(tool):
            raise SystemExit(f"{tool} is not on PATH")
    check_same_filesystem(args.input_dir, args.state_dir)
    sra_ids = sorted(subdirs(os.path.join(args.state_dir, "downloaded")))
    validations_csv = os.path.join(args.state_dir, "validations.csv")
    print(f"{len(sra_ids)} downloads to check with {args.threads} threads", flush=True)

    results, ok_bytes = Counter(), 0
    with ThreadPoolExecutor(args.threads) as pool:
        futures = [pool.submit(check_download, sra_id, args) for sra_id in sra_ids]
        for future in as_completed(futures):
            record = future.result()
            append_record(validations_csv, VALIDATION_COLUMNS, record)
            results[record["result"]] += 1
            ok_bytes += record["bytes"] if record["result"] == "ok" else 0
            print(f"{record['sra_id']}  {record['result']}  {record['seconds']} s  "
                  f"{record['bytes'] / 1e6:.1f} MB  {record['reason']}", flush=True)
    print(f"checked {len(sra_ids)}: {results['ok']} ok ({ok_bytes / 1e9:.1f} GB moved to the input dir), "
          f"{results['invalid']} invalid", flush=True)


# --- status ---

def print_method_row(name, rows):
    results = Counter(row["result"] for row in rows)
    ok = [row for row in rows if row["result"] == "ok"]
    ok_bytes = sum(int(row["bytes"]) for row in ok)
    ok_seconds = sum(int(row["seconds"]) for row in ok)
    speed = f"{ok_bytes / 1e6 / ok_seconds:.1f}" if ok_seconds else "-"
    counts = "".join(f" {results[result]:>{max(len(result), 6)}}" for result in ATTEMPT_RESULTS)
    print(f"{name:<12} {len(rows):>8}{counts} {ok_bytes / 1e9:>8.1f} {speed:>6}")


def print_failure_reasons(title, rows):
    """Print the 10 most frequent reasons of the counted failures.

    The run accession becomes <run>, so that the same cause groups.
    """
    counted = [row for row in rows if row["result"] in COUNTED_RESULTS]
    reasons = Counter(RUN_ACCESSION.sub("<run>", row["reason"]) for row in counted)
    print(f"\n{title}: {sum(reasons.values())}")
    for reason, count in reasons.most_common(10):
        print(f"{count:>7}  {reason}")


def status(args):
    attempt_rows = read_records(os.path.join(args.state_dir, "attempts.csv"))
    validation_rows = read_records(os.path.join(args.state_dir, "validations.csv"))
    if args.job:
        attempt_rows = [row for row in attempt_rows if row["job_id"] == args.job]
        validation_rows = [row for row in validation_rows if row["job_id"] == args.job]

    n_samples = len({row["sra_id"] for row in attempt_rows})
    n_jobs = len({row["job_id"] for row in attempt_rows})
    print(f"attempts.csv: attempts {len(attempt_rows)}, samples {n_samples}, job IDs {n_jobs}\n")
    header = "".join(f" {result:>{max(len(result), 6)}}" for result in ATTEMPT_RESULTS)
    print(f"{'method':<12} {'attempts':>8}{header} {'ok GB':>8} {'MB/s':>6}")
    for method in sorted({row["method"] for row in attempt_rows}):
        print_method_row(method, [row for row in attempt_rows if row["method"] == method])
    print_method_row("total", attempt_rows)

    ok_checks = [row for row in validation_rows if row["result"] == "ok"]
    moved_gb = sum(int(row["bytes"]) for row in ok_checks) / 1e9
    print(f"\nvalidations.csv: {len(validation_rows)} content checks, {len(ok_checks)} ok "
          f"({moved_gb:.1f} GB moved to the input dir), {len(validation_rows) - len(ok_checks)} invalid")

    print_failure_reasons("failure reasons of the attempts", attempt_rows)
    print_failure_reasons("failure reasons of the content checks", validation_rows)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["runinfo", "plan", "run", "validate", "status"])
    parser.add_argument("--ids", help="runinfo, plan, run: CSV with an sra_id column and an optional Bytes column")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--state-dir", help="default: <input-dir>.download")
    parser.add_argument("--max-attempts", type=int, default=3)
    parser.add_argument("--methods", default="aws-http:4,ena-ftp:4",
                        help="run: kingfisher methods and their number of threads")
    parser.add_argument("--attempt-timeout", type=int, default=6 * 3600, help="run: seconds before an attempt is stopped")
    parser.add_argument("--threads", type=int, default=4, help="validate: downloads checked at the same time")
    parser.add_argument("--file-timeout", type=int, default=2 * 3600, help="validate: seconds for the check of one file")
    parser.add_argument("--job", help="status: only the rows of this job ID (a download job or a validate job)")
    args = parser.parse_args()
    if args.command in ("runinfo", "plan", "run") and not args.ids:
        parser.error(f"{args.command} needs --ids")
    args.state_dir = args.state_dir or os.path.normpath(args.input_dir) + ".download"
    args.methods = parse_methods(args.methods)
    return args


def main():
    args = parse_args()
    if args.command == "runinfo":
        fetch_runinfo(args)
    elif args.command == "plan":
        _, failures = read_history(args.state_dir)
        plan(args, failures)
    elif args.command == "run":
        run(args)
    elif args.command == "validate":
        check_downloads(args)
    elif args.command == "status":
        status(args)


if __name__ == "__main__":
    main()

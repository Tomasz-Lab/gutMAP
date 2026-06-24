#!/bin/bash -l
#SBATCH -J kingfisher
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=50GB
#SBATCH --time=48:00:00
#SBATCH -A plggutmap100k-cpu
#SBATCH -p plgrid
#SBATCH --output="/net/scratch/hscra/plgrid/plgpkica/metagenome_proj/results/kf/output/%J.out"
#SBATCH --error="/net/scratch/hscra/plgrid/plgpkica/metagenome_proj/results/kf/error/%J.err"
#
# Single-job kingfisher downloader.
#
# The bash part is just SLURM directives + conda env activation. All queue
# logic and worker orchestration is Python: an in-memory queue.Queue is shared
# across N_AWS + N_PREFETCH + N_ENA daemon threads, each blocking on its own
# kingfisher subprocess. queue.Queue.get_nowait() is thread-safe under the
# GIL, so no flock/lockfile is needed. Each worker class differs only in
# which kingfisher -m method is tried FIRST (kingfisher falls back through
# the rest on failure), spreading load across AWS, NCBI prefetch and ENA-FTP.


set -o pipefail

# ============================================================
# Worker mix - tune these without touching anything else.
# Total concurrent kingfishers = N_AWS + N_PREFETCH + N_ENA.
# ============================================================
export N_AWS=${N_AWS:-3}            # method order: aws-http  prefetch ena-ftp
export N_PREFETCH=${N_PREFETCH:-2}  # method order: prefetch  aws-http ena-ftp
export N_ENA=${N_ENA:-3}            # method order: ena-ftp   prefetch aws-http
export KF_THREADS=${KF_THREADS:-1}


# ============================================================
# Paths
# ============================================================
PROJ=/net/scratch/hscra/plgrid/plgpkica/metagenome_proj
export MASTER_CSV="$PROJ/data/input_dataset.csv"
export INPUT_DIR="$PLG_GROUPS_STORAGE/plgggutmap/input_mg/input_dataset"
export PRESENCE_DIR="$INPUT_DIR"   # checked for "already downloaded" sra_ids
export GUTMAP_BASE="$PLG_GROUPS_STORAGE/plgggutmap/output_mg/gutmap"
export CANONICAL_STEP="02_fastp_qc"
export SKIP_STATES="valid"          # space-separated, e.g. "valid invalid"

export RUN_DIR="$PROJ/results/kf/run_${SLURM_JOB_ID:-local}"
mkdir -p "$RUN_DIR" "$INPUT_DIR"

module load Miniconda3/23.9.0-0
source activate /net/scratch/hscra/plgrid/plgpkica/metagenome_proj/conda/env-53759597ff597e1f50216b1bb7cc4db4

#export NCBI_API_KEY=""

# ============================================================
# Everything else: in-process Python orchestrator
# ============================================================
python <<'PY'
import csv
import os
import queue
import subprocess
import sys
import threading
from pathlib import Path

MASTER       = os.environ["MASTER_CSV"]
PRESENCE_DIR = os.environ["PRESENCE_DIR"]
INPUT_DIR    = os.environ["INPUT_DIR"]
RUN_DIR      = Path(os.environ["RUN_DIR"])
N_AWS        = int(os.environ["N_AWS"])
N_PREFETCH   = int(os.environ["N_PREFETCH"])
N_ENA        = int(os.environ["N_ENA"])
KF_THREADS   = int(os.environ["KF_THREADS"])

# ---- build queue ---------------------------------------------------------
done = set(os.listdir(PRESENCE_DIR)) if os.path.isdir(PRESENCE_DIR) else set()
on_disk = len(done)

# Treat already-classified samples (gutmap/<state>/<canonical_step>/<sra_id>)
# as done even if their input files were deleted to free space. Without this,
# removing a valid sample's input would re-queue it for download.
GUTMAP_BASE   = os.environ.get("GUTMAP_BASE", "")
CANONICAL_STEP = os.environ.get("CANONICAL_STEP", "02_fastp_qc")
SKIP_STATES   = os.environ.get("SKIP_STATES", "").split()
classified = set()
for st in SKIP_STATES:
    d = os.path.join(GUTMAP_BASE, st, CANONICAL_STEP)
    if os.path.isdir(d):
        classified |= set(os.listdir(d))
done |= classified

remaining = []
with open(MASTER) as f:
    for row in csv.DictReader(f):
        sid = (row.get("sra_id") or "").strip()
        if sid and sid not in done:
            remaining.append(sid)

print(f"[work] on_disk={on_disk}  classified_skip={len(classified)}  "
      f"({' '.join(SKIP_STATES) or 'none'})  remaining={len(remaining)}", flush=True)
if not remaining:
    print("Nothing to download. Exiting.")
    sys.exit(0)

q = queue.Queue()
for sid in remaining:
    q.put(sid)

# ---- worker --------------------------------------------------------------
def worker(tag: str, methods: list[str]) -> None:
    log_path = RUN_DIR / f"worker_{tag}.log"
    with open(log_path, "w", buffering=1) as logf:   # line-buffered
        logf.write(f"[{tag}] start (methods: {' '.join(methods)})\n")
        while True:
            try:
                sra_id = q.get_nowait()
            except queue.Empty:
                logf.write(f"[{tag}] queue empty, exiting\n")
                return

            dest = os.path.join(INPUT_DIR, sra_id)
            if os.path.isdir(dest) and os.listdir(dest):
                logf.write(f"[{tag}] skip {sra_id} (already present)\n")
                continue

            logf.write(f"[{tag}] start {sra_id}\n")
            cmd = [
                "kingfisher", "get",
                "-r", sra_id,
                "-t", str(KF_THREADS),
                "-f", "fastq.gz", "fastq", "sra",
                "-m", *methods,
                "--output-directory", dest,
                "--check-md5sums",
                "--guess-aws-location",
                "--debug",   # kingfisher's verbose native logging -> stderr -> worker log
            ]
            # PYTHONUNBUFFERED=1 forces line-buffered stderr/stdout in the child
            # so kingfisher's --debug output streams live into the worker log
            # instead of sitting in its ~4 KB block buffer until exit.
            env = {**os.environ, "PYTHONUNBUFFERED": "1"}
            r = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT, env=env)
            if r.returncode == 0:
                logf.write(f"[{tag}] OK   {sra_id}\n")
            else:
                logf.write(f"[{tag}] FAIL {sra_id} (rc={r.returncode})\n")

# ---- spawn ---------------------------------------------------------------
threads: list[threading.Thread] = []

def spawn(n: int, tag: str, methods: list[str]) -> None:
    for i in range(1, n + 1):
        t = threading.Thread(target=worker, args=(f"{tag}-{i}", methods))
        t.start()
        threads.append(t)

spawn(N_AWS,      "aws",      ["aws-http", "prefetch", "ena-ftp"])
spawn(N_PREFETCH, "prefetch", ["prefetch", "aws-http", "ena-ftp"])
spawn(N_ENA,      "ena",      ["ena-ftp",  "prefetch", "aws-http"])

print(f"Spawned {len(threads)} workers (aws={N_AWS}, prefetch={N_PREFETCH}, ena={N_ENA})")
print(f"Worker logs: {RUN_DIR}/worker_*.log")

# ---- wait + summarize ----------------------------------------------------
for t in threads:
    t.join()

print()
print("All workers complete.")
final = len(os.listdir(INPUT_DIR)) if os.path.isdir(INPUT_DIR) else 0
print(f"Final on-disk: {final} subdirs in {INPUT_DIR}")
print(f"Queue pending: {q.qsize()}")
PY

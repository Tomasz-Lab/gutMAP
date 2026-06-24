#!/bin/bash -l
#SBATCH -J workdir_size
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=4GB
#SBATCH --time=1:00:00
#SBATCH -A plggutmap100k-cpu
#SBATCH -p plgrid
#SBATCH --output="/net/scratch/hscra/plgrid/plgpkica/metagenome_proj/logs/measure/%J.out"
#SBATCH --error="/net/scratch/hscra/plgrid/plgpkica/metagenome_proj/logs/measure/%J.err"
#
# Per-batch workdir sizing as a single SLURM job.
#
# Usage (invoked by run_batches.sh with `sbatch --wait`):
#   sbatch --wait measure_workdirs.sbatch <nf_log> <batch_name>
#
# 1. Parses the Nextflow log for `submitted process ... workDir: ...` events.
# 2. Runs `du -sb --apparent-size` on each unique workdir in parallel
#    (ProcessPoolExecutor, 10 workers = --cpus-per-task).
# 3. Appends rows to the shared CSV at
#    $PROJ/logs/workdir_sizes_batches.csv
#    with columns: batch,sra_id,pipeline_step,workdir_path,workdir_size
#
# No append-lock is needed: run_batches.sh uses sbatch --wait so at most one
# measurement job is active at a time.

set -euo pipefail

NF_LOG="${1:?usage: measure_workdirs.sbatch <nf_log> <batch_name>}"
BATCH_NAME="${2:?usage: measure_workdirs.sbatch <nf_log> <batch_name>}"

PROJ=/net/scratch/hscra/plgrid/plgpkica/metagenome_proj
SHARED_CSV="$PROJ/logs/workdir_sizes_batches.csv"
mkdir -p "$(dirname "$SHARED_CSV")" "$PROJ/logs/measure"

export NF_LOG BATCH_NAME SHARED_CSV

python3 <<'PY'
import csv
import os
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

NF_LOG     = Path(os.environ["NF_LOG"])
BATCH_NAME = os.environ["BATCH_NAME"]
SHARED_CSV = Path(os.environ["SHARED_CSV"])

# --- 1. parse log (same regex as the old extract_workdirs.py) ---
LINE_RE = re.compile(
    r"submitted process (?P<step>\S+) \((?P<tag>[^)]*)\) > jobId: \d+; workDir: (?P<workdir>\S+)"
)
SRA_RE = re.compile(r"[A-Z]{3}\d+")

rows = []
with NF_LOG.open(encoding="utf-8", errors="replace") as f:
    for line in f:
        m = LINE_RE.search(line)
        if not m:
            continue
        step, tag, workdir = m["step"], m["tag"], m["workdir"]
        if step == "BAKTA_BATCH":
            sra_id = ""   # filled in during sizing, one row per `.bakta/` member
        else:
            ids = SRA_RE.findall(tag)
            sra_id = ids[-1] if ids else ""
        rows.append({
            "sra_id":        sra_id,
            "pipeline_step": step,
            "workdir_path":  workdir,
        })

bakta_workdirs    = sorted({r["workdir_path"] for r in rows if r["pipeline_step"] == "BAKTA_BATCH"})
nonbakta_workdirs = sorted({r["workdir_path"] for r in rows if r["pipeline_step"] != "BAKTA_BATCH"})
print(f"parsed {len(rows)} rows, {len(nonbakta_workdirs)} non-bakta + {len(bakta_workdirs)} bakta workdirs from {NF_LOG.name}", flush=True)


def du_path(p: str) -> str:
    """Returns size in bytes as a string, or 'MISSING'/'ERROR'."""
    try:
        r = subprocess.run(
            ["du", "-sb", "--apparent-size", p],
            capture_output=True, text=True, timeout=900,
        )
        if r.returncode != 0:
            return "MISSING"
        return r.stdout.split("\t", 1)[0]
    except subprocess.TimeoutExpired:
        return "ERROR"
    except Exception:
        return "ERROR"


def size_workdir(args):
    """Returns (workdir, [(sra_id_or_None, size_str), ...]).
    - Non-bakta: list has exactly one entry; sra_id is None (caller uses upstream id).
    - Bakta: one entry per `<sra>.bakta/` subdir found, size = du(.bakta) + size(.stdout.log).
    """
    workdir, is_bakta = args
    if not is_bakta:
        return (workdir, [(None, du_path(workdir))])
    # Bakta: expand
    if not os.path.isdir(workdir):
        return (workdir, [("", "MISSING")])
    out = []
    for e in sorted(os.listdir(workdir)):
        if not e.endswith(".bakta"):
            continue
        member_dir = os.path.join(workdir, e)
        if not os.path.isdir(member_dir):
            continue
        sra = e[:-len(".bakta")]
        log_path = os.path.join(workdir, f"{sra}.bakta.stdout.log")
        dsize = du_path(member_dir)
        if dsize in ("MISSING", "ERROR"):
            out.append((sra, dsize))
            continue
        try:
            lsize = os.path.getsize(log_path) if os.path.isfile(log_path) else 0
        except OSError:
            lsize = 0
        out.append((sra, str(int(dsize) + lsize)))
    if not out:
        out.append(("", "MISSING"))
    return (workdir, out)


# --- 2. measure in parallel (10 workers = --cpus-per-task) ---
tasks = [(w, False) for w in nonbakta_workdirs] + [(w, True) for w in bakta_workdirs]
results: dict[str, list[tuple[str | None, str]]] = {}
with ProcessPoolExecutor(max_workers=10) as pool:
    for workdir, members in pool.map(size_workdir, tasks):
        results[workdir] = members

# --- 3. assemble output rows ---
out_rows = []
for r in rows:
    wd = r["workdir_path"]
    members = results.get(wd, [(None, "MISSING")])
    if r["pipeline_step"] == "BAKTA_BATCH":
        # Expand: one upstream row -> N member rows (one per .bakta/ dir)
        for sra, size in members:
            out_rows.append({
                "batch":         BATCH_NAME,
                "sra_id":        sra or "",
                "pipeline_step": "BAKTA_BATCH",
                "workdir_path":  wd,
                "workdir_size":  size,
            })
    else:
        size = members[0][1]
        out_rows.append({
            "batch":         BATCH_NAME,
            "sra_id":        r["sra_id"],
            "pipeline_step": r["pipeline_step"],
            "workdir_path":  wd,
            "workdir_size":  size,
        })

# --- 4. append to shared CSV (write header iff file is new/empty) ---
fields = ["batch", "sra_id", "pipeline_step", "workdir_path", "workdir_size"]
new_file = not SHARED_CSV.exists() or SHARED_CSV.stat().st_size == 0
with SHARED_CSV.open("a", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    if new_file:
        w.writeheader()
    w.writerows(out_rows)

n_missing = sum(1 for r in out_rows if r["workdir_size"] == "MISSING")
n_error   = sum(1 for r in out_rows if r["workdir_size"] == "ERROR")
n_ok      = len(out_rows) - n_missing - n_error
print(f"appended {len(out_rows)} rows -> {SHARED_CSV}", flush=True)
print(f"  OK      : {n_ok}", flush=True)
print(f"  MISSING : {n_missing}", flush=True)
print(f"  ERROR   : {n_error}", flush=True)
PY

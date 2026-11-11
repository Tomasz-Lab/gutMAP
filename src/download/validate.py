"""Check that a sample dir holds a complete download that the pipeline can read.

The dir name is the sra_id. The file-set check only lists the files. The content
check also reads each file with vdb-validate or gzip -t.
"""
import os
import subprocess


def pipeline_file_sets(sra_id):
    """The file sets that mtg_pipeline.nf reads (the input step at the top of the workflow)."""
    sra = f"{sra_id}.sra"
    r1, r2, single = (f"{sra_id}{suffix}.fastq.gz" for suffix in ("_1", "_2", ""))
    return [{sra}, {r1, r2}, {r1, r2, single}, {single}]


def file_set_problems(sample_dir):
    sra_id = os.path.basename(os.path.normpath(sample_dir))
    entries = list(os.scandir(sample_dir))
    names = {e.name for e in entries}
    problems = []
    if names not in pipeline_file_sets(sra_id):
        problems.append("files not readable by the pipeline: " + (" ".join(sorted(names)) or "none"))
    for e in entries:
        if not e.is_file(follow_symlinks=False):
            problems.append(f"not a regular file: {e.name}")
        elif e.stat().st_size == 0:
            problems.append(f"empty file: {e.name}")
    return problems


def content_problem(path, timeout):
    """Return None if the reader accepts the file, else the reason."""
    cmd = ["vdb-validate", path] if path.endswith(".sra") else ["gzip", "-t", path]
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           text=True, errors="replace", timeout=timeout)
    except subprocess.TimeoutExpired:
        return f"{os.path.basename(path)}: {cmd[0]} timeout after {timeout} s"
    if p.returncode == 0:
        # vdb-validate also exits with 0 for a file that is not an SRA archive, and then prints nothing.
        if path.endswith(".sra") and "is consistent" not in p.stdout:
            return f"{os.path.basename(path)}: vdb-validate found no SRA archive"
        return None
    lines = [line.strip() for line in p.stdout.splitlines() if line.strip()]
    error = next((line for line in lines if "err" in line.lower()), lines[-1] if lines else "")
    error = error.replace(f"{path}: ", "").replace(path, os.path.basename(path))
    return f"{os.path.basename(path)}: {cmd[0]} exit {p.returncode}: {error[:200]}"


def content_problems(sample_dir, timeout):
    paths = sorted(os.path.join(sample_dir, name) for name in os.listdir(sample_dir))
    return [problem for problem in (content_problem(path, timeout) for path in paths) if problem]


def find_problems(sample_dir, content=False, timeout=7200):
    """Return a list of problems. An empty list means that the dir passes."""
    problems = file_set_problems(sample_dir)
    if content and not problems:
        problems = content_problems(sample_dir, timeout)
    return problems

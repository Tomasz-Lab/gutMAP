#!/usr/bin/env python3
"""Collect the Slurm job history of Helios into one CSV. Run it on Helios.

The script calls `hpc-jobs-history --verbose` for one day at a time, from 00:00
to 00:00 of the next day, and writes one row for each job. hpc-jobs-history
selects jobs by their end time (Helios local time). When the CSV has
rows, the script starts on the day before the newest end time in it. A job that
is already in the CSV gets the new row.
If a call fails, the script keeps the rows of the days before it. Run it again
to continue.

    python3 collect_hjh.py --out hpc-jobs-history.csv --start 2026-03-01
    python3 collect_hjh.py --out hpc-jobs-history.csv

Python 3.7 or later, standard library only.
"""

import argparse
import csv
import io
import os
import re
import subprocess
import sys
from datetime import datetime, timedelta, timezone

FIELDS = [
    "job_id", "partition", "Account", "pipeline_step", "sra_id", "nodes", "cores",
    "decl_mem[B]", "Max._node_mem.", "mem_usage[%]", "efficiency[%]", "cpu_used[s]",
    "billing[CPUh]", "wall_used[s]", "wall_req[s]", "In_queue", "Exit", "State", "end_time",
]
TIME_FORMAT = "%Y-%m-%d %H:%M:%S"
MEM_UNITS = {
    "B": 1, "KiB": 2**10, "MiB": 2**20, "GiB": 2**30, "TiB": 2**40,
    "KB": 10**3, "MB": 10**6, "GB": 10**9, "TB": 10**12,
}


def is_missing(value):
    return re.fullmatch(r"-*", value) is not None


def to_seconds(value):
    """'[D-]HH:MM:SS' -> seconds."""
    if is_missing(value):
        return ""
    match = re.fullmatch(r"(?:(\d+)-)?(\d+):(\d+):(\d+)", value)
    if not match:
        return value
    days, hours, minutes, seconds = (int(x or 0) for x in match.groups())
    return str(((days * 24 + hours) * 60 + minutes) * 60 + seconds)


def to_bytes(value):
    """'8.0GiB' -> bytes."""
    match = re.fullmatch(r"([0-9.]+)\s*([A-Za-z]+)", value)
    if not match or match.group(2) not in MEM_UNITS:
        return value
    return str(int(float(match.group(1)) * MEM_UNITS[match.group(2)]))


def to_percent(value):
    value = value.rstrip("%")
    return "" if is_missing(value) else value


def to_cpuh(value):
    """'0.55 CPUh' -> '0.55'."""
    value = value.replace("CPUh", "").strip()
    return "" if is_missing(value) else value


def to_unix_time(value):
    """Helios local time, stored as if it were UTC, as in the CSVs from before."""
    try:
        moment = datetime.strptime(value, TIME_FORMAT)
    except ValueError:
        return value
    return str(int(moment.replace(tzinfo=timezone.utc).timestamp()))


def split_name(name):
    """'nf-FASTP_QC_(fastp__ERR866589)' -> ('FASTP_QC', 'ERR866589')."""
    step_match = re.match(r"nf-(.+?)_\(", name)
    sra_match = re.search(r"([A-Z]{3}\d+)\)?$", name)
    if step_match:
        step = step_match.group(1)
    elif name.startswith("nf-"):
        step = name[len("nf-"):]
    else:
        step = name
    return step, sra_match.group(1) if sra_match else ""


def to_row(job):
    """One job of hpc-jobs-history -> the values of FIELDS."""
    step, sra_id = split_name(job["Name"])
    wall_used = to_seconds(job["Wall._Used"])
    if "CPU._used" in job:
        cpu_used = to_seconds(job["CPU._used"])
    elif job["Cores"].isdigit() and wall_used.isdigit():
        cpu_used = str(int(job["Cores"]) * int(wall_used))
    else:
        cpu_used = ""
    return [
        job["ID"], job["Partition"], job["Account"], step, sra_id, job["Nodes"], job["Cores"],
        to_bytes(job["Decl._mem"]), job["Max._node_mem."], to_percent(job["Mem._%_usage"]),
        to_percent(job["Eff."]), cpu_used, to_cpuh(job.get("Billing", "")), wall_used,
        to_seconds(job["Wall._Req."]), job["In_queue"], job["Exit"], job["State"],
        to_unix_time(job["End_Time"]),
    ]


def parse_table(text):
    """The table of hpc-jobs-history --verbose -> (jobs as dicts, number of bad lines).

    Columns are separated by 2 or more spaces. In the header, a single space
    becomes '_' ('Decl. mem' -> 'Decl._mem').
    """
    header, jobs, bad_lines = None, [], 0
    for line in text.splitlines():
        if not line.strip() or re.fullmatch(r"[\s-]+", line):
            continue
        columns = re.split(r"\s{2,}", line.strip())
        if header is None:
            header = [column.replace(" ", "_") for column in columns]
        elif len(columns) == len(header):
            jobs.append(dict(zip(header, columns)))
        else:
            bad_lines += 1
    return jobs, bad_lines


def fetch_day(command, start, end):
    """hpc-jobs-history for [start, end) -> {job_id: CSV line}."""
    output = subprocess.run(
        [command, "--verbose", "--start", start.strftime(TIME_FORMAT), "--end", end.strftime(TIME_FORMAT)],
        check=True, stdout=subprocess.PIPE, text=True,
    ).stdout
    jobs, bad_lines = parse_table(output)
    lines = {}
    for job in jobs:
        try:
            row = to_row(job)
        except KeyError as error:
            sys.exit(f"hpc-jobs-history has no column {error}. Its header: {list(job)}")
        buffer = io.StringIO()
        csv.writer(buffer, lineterminator="\n").writerow(row)
        lines[row[0]] = buffer.getvalue()
    return lines, bad_lines


def newest_end_time(path):
    """The newest end_time in the CSV, or None if it has no rows."""
    if not os.path.exists(path):
        return None
    with open(path, newline="") as f:
        reader = csv.reader(f)
        if next(reader, None) != FIELDS:
            sys.exit(f"{path} does not have the columns of this script")
        newest = max((int(row[-1]) for row in reader if row[-1].isdigit()), default=None)
    if newest is None:
        return None
    return datetime.fromtimestamp(newest, timezone.utc).replace(tzinfo=None)


def write_csv(path, new_lines):
    """Old rows of jobs that are not in new_lines, then new_lines. Atomic."""
    temporary = path + ".tmp"
    with open(temporary, "w", newline="") as out:
        out.write(",".join(FIELDS) + "\n")
        if os.path.exists(path):
            with open(path, newline="") as f:
                reader = csv.reader(f)
                next(reader)
                writer = csv.writer(out, lineterminator="\n")
                for row in reader:
                    if row[0] not in new_lines:
                        writer.writerow(row)
        out.writelines(new_lines.values())
    os.replace(temporary, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", required=True, help="the CSV to update")
    parser.add_argument("--start", help="YYYY-MM-DD, used only when the CSV has no rows")
    parser.add_argument("--command", default="hpc-jobs-history", help="for tests with a stub")
    args = parser.parse_args()

    newest = newest_end_time(args.out)
    if newest:
        first_day = newest.date() - timedelta(days=1)
    elif args.start:
        first_day = datetime.strptime(args.start, "%Y-%m-%d").date()
    else:
        sys.exit(f"{args.out} has no rows: give --start")
    today = datetime.now().date()

    new_lines, day, failed = {}, first_day, False
    print(f"collect from {first_day} to {today}")
    while day <= today:
        midnight = datetime.combine(day, datetime.min.time())
        try:
            lines, bad_lines = fetch_day(args.command, midnight, midnight + timedelta(days=1))
        except subprocess.CalledProcessError:
            failed = True
            break
        new_lines.update(lines)
        print(f"{day}: {len(lines)} jobs" + (f", {bad_lines} bad lines" if bad_lines else ""))
        day += timedelta(days=1)

    write_csv(args.out, new_lines)
    print(f"{len(new_lines)} jobs from {args.command} -> {args.out}")
    if failed:
        sys.exit(f"{args.command} failed for {day}. Run the script again to continue.")


if __name__ == "__main__":
    main()

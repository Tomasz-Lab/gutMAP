# Download of input data

`download.py` downloads SRA and ENA runs with kingfisher into the input dir of the pipeline. Only downloads that pass all checks go into the input dir.

## Files

- `download.py`: the commands `runinfo`, `plan`, `run`, `validate` and `status`. For the options, run `python3 download.py --help`.
- `runinfo.py`: the requests to NCBI E-utilities for `runinfo`.
- `validate.py`: the checks of one sample dir. `run` and `validate` use the same checks.
- `download.sbatch`, `validate.sbatch`: the Slurm jobs for `run` and `validate`.

## Data flow

```
<input-dir>.download/staging/<id>      kingfisher output during an attempt
<input-dir>.download/downloaded/<id>   correct file set, waits for the content check
<input-dir>/<id>                       content check ok, the pipeline reads it
```

- File-set check (`run`): the files are `<id>.sra`, or `<id>_1.fastq.gz` and `<id>_2.fastq.gz` (with or without `<id>.fastq.gz`), or `<id>.fastq.gz`. No file is empty.
- Content check (`validate`): `vdb-validate` for `.sra` and `gzip -t` for `.fastq.gz`.
- A sample with `--max-attempts` failures (default 3) is exhausted. `run` does not try it again.

The state dir `<input-dir>.download/` also holds `attempts.csv` (one row for each download attempt), `validations.csv` (one row for each content check) and `logs/` (the kingfisher log of each failed attempt).

## Runinfo

`runinfo` adds the NCBI SRA runinfo row of each new sample to `<input-dir>.download/runinfo.csv`. The batch selection uses `bases` and `spots` from it to estimate the disk use of each sample. `plan` uses `size_MB` when the `--ids` CSV has no `Bytes` value.

- Run it before `plan`, on a login node. It sends 200 IDs in each request, so 100k samples take some minutes.
- It asks only for the samples that are not in `runinfo.csv`. If it stops, run it again.
- If a request for 200 samples still fails after its retries, `runinfo` asks for each sample of it alone. NCBI fails on some samples every time, for example `SRR769534`.
- The samples without runinfo go to `runinfo_missing.csv`, with the reason: not returned, or the error.
- Without `NCBI_API_KEY`, NCBI allows 3 requests each second. With a key, 10. Set the key as an environment variable. `NCBI_EMAIL` and `NCBI_TOOL` are optional.
- efetch sends no header row. If a row does not have the expected 47 fields, `runinfo` stops.

## Requirements

- Python 3.7 or later.
- A conda env with kingfisher and sra-tools: `conda create -p ENV -c conda-forge -c bioconda kingfisher=0.5.0 sra-tools=3.4.1`. Set `KF_ENV` to `ENV`.

## Run with Slurm

The `Makefile` at the repo root runs these commands with the settings of `config.mk`. Run `make` from the repo root:

1. Set `IDS`, `INPUT_DIR`, `KF_ENV`, `ACCOUNT` and `PARTITION` in `config.mk`.
2. Get the runinfo, show the plan, submit the download job and the validate job, and show the result:

   ```
   make runinfo
   make download-plan
   make download
   make download-status
   ```

- `make download` stops if a batches job is in the queue. Downloads and batches must not run at the same time.
- `make validate` submits only the validate job.
- `make help` lists all targets. A setting on the command line replaces the one in `config.mk`, for example `make download METHODS=ena-ftp:4`.

Without `make`, submit the jobs from this dir:

```
python3 download.py runinfo --ids samples.csv --input-dir INPUT_DIR
python3 download.py plan --ids samples.csv --input-dir INPUT_DIR
export KF_ENV=ENV SBATCH_ACCOUNT=ACCOUNT SBATCH_PARTITION=PARTITION
J=$(sbatch --parsable download.sbatch --ids samples.csv --input-dir INPUT_DIR)
sbatch --dependency=singleton,afterany:$J validate.sbatch --input-dir INPUT_DIR
python3 download.py status --input-dir INPUT_DIR
```

With `--dependency=singleton`, only one download job and one validate job run at a time.

## Run without Slurm

```
export PATH=$KF_ENV/bin:$PATH
python3 download.py run --ids samples.csv --input-dir INPUT_DIR > download.log 2>&1
python3 download.py validate --input-dir INPUT_DIR
```

Start only one `run` at a time for each input dir. A second `run` removes the `staging/` dir of the first.

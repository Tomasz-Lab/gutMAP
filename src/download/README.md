# Download of input data

`download.py` downloads SRA and ENA runs with kingfisher into the input dir of the pipeline. Only downloads that pass all checks go into the input dir.

## Files

- `download.py`: the commands `plan`, `run`, `validate` and `status`. For the options, run `python3 download.py --help`.
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

## Requirements

- Python 3.7 or later.
- A conda env with kingfisher and sra-tools: `conda create -p ENV -c conda-forge -c bioconda kingfisher=0.5.0 sra-tools=3.4.1`. Set `KF_ENV` to `ENV`.

## Run with Slurm

Submit the jobs from this dir:

```
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

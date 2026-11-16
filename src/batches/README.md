# Batch runs of the pipeline

`batches.py` runs `mtg_pipeline.nf` in batches that fit a disk budget. It records the state of each sample in one table and moves the accepted results to the final dataset.

## Files

- `batches.py`: the commands `plan`, `run`, `status`, `sample` and `promote`. For the options, run `python3 batches.py <command> --help`.
- `batches.sbatch`: the Slurm job for `run`. Nextflow runs in this job and submits the tasks.

## Data flow

```
<input-dir>/<id>/                      downloaded input (from ../download)
<input-dir>.download/runinfo.csv       NCBI runinfo, gives the size of each sample

<staging-dir>/samples.tsv              the state of each sample: the only record
<staging-dir>/traces.tsv               the Nextflow traces of all batches, one table
<staging-dir>/promote.csv              one row for each dir that promote moved
<staging-dir>/not_planned.tsv          the samples in state too large or no size, from the last plan or run
<staging-dir>/batch_NNNN/
    batch.csv                          the samples of the batch and their disk reserve
    nextflow.log
    pipeline_info/<run>/execution_trace.tsv
    published/<step>/<id>/             the outputs of the pipeline, until promote; only the logs
                                       of the samples that promote does not take
    failed_logs/<process>/<id>/        the logs of the failed tasks

<valid-dir>/<step>/<id>/               accepted results, without logs
<logs-dir>/<step>/<id>/logs/           the logs of the accepted samples
```

The input dir, the staging dir, the valid dir and the logs dir must be on one filesystem, because `promote` moves the dirs with `os.rename`.

## Sample states

| State | Meaning | Next action |
|---|---|---|
| `eligible` | Downloaded and not in a batch yet, or `stopped` | `run` puts it in a batch |
| `not downloaded` | Not in the input dir | Download it |
| `too large` | Downloaded, but its disk reserve is larger than `--disk-budget` | Increase the budget, or leave it out |
| `no size` | Downloaded, but runinfo has no bases for it and its input dir is empty | Download it again |
| `running` | In the batch that runs now, or in a batch whose run stopped | The next `run` finishes the batch |
| `done` | No task failed and the size gate passed | `promote` |
| `failed` | A task failed, or the size gate failed. This is final. | `promote` if the gate passed, else none |
| `stopped` | The run stopped before the sample ended (reason `run stopped`) | `run` puts it in a later batch |
| `promoted` | The results are in the valid dir | None |

The size gate is the acceptance rule of the dataset: the sylph profile is larger than 218 bytes (more than the header) and the contigs file is not empty. A step that did not run because its input was empty, for example no bins, is not a failure.

Only Nextflow retries a task, in the same run. A process that needs a retry sets it in its own code, with more resources for the retry, for example `FILTER_AND_ANNOTATE`. All other failed tasks are final (`errorStrategy = 'ignore'`). A sample does not go into a later batch because a task failed.

If the run stopped, a sample gets `stopped` when one of its tasks has the status `ABORTED`, or when no task failed and an endpoint (`SYLPH_TAXONOMY`, `CALCULATE_TPM_AND_ANNOTATE`, `FILTER_AND_ANNOTATE`) did not run. All other samples of the batch get `done` or `failed`, as after a normal run. A sample whose branch a filter skipped (empty reads, no contigs) also gets `stopped`, and runs again.

After each batch, `run` removes from `published/` the outputs of the samples that `promote` does not take: `stopped`, and gate `no`. The logs stay, and so do the dirs of `02_fastp_qc` and `03_non_human_reads`.

`samples.tsv` also has the sizes of the sylph profile and the contigs (`sylph_bytes`, `contigs_bytes`), and one column for each pipeline step: `ok`, `failed <exit code>`, or `-` when the step did not run.

## Disk budget

Each sample reserves `--reserve-factor` (default 4.5) times its FASTQ size. The FASTQ size comes from runinfo: 2 × `bases` + 120 B × `spots`. If runinfo has no bases for a sample, the sample reserves `--input-factor` (default 36) times the size of its input dir. A batch takes the eligible samples in the order of `--ids` until the next sample does not fit in `--disk-budget`.

A sample whose reserve is larger than the budget gets the state `too large`, and no batch takes it. `plan`, `run` and `status` calculate this state again at each call from the budget that you give. `plan` and `run` write these samples and their reserve to `not_planned.tsv`. When no sample is `eligible`, `run` stops.

Downloads must not run at the same time as a batch, because the budget does not include them.

## Run on Helios

The `Makefile` at the repo root runs these commands with the settings of `config.mk`. See the README at the repo root.

Run the commands from the launch dir, the dir that holds `src/` and the launch-dir `nextflow.config`.

1. Download the samples and their runinfo with `../download`.
2. Show the next batch:

   ```
   python3 src/batches/batches.py plan --staging-dir S --ids IDS.csv --input-dir INPUT_DIR --disk-budget 20T
   ```

3. Submit the batches:

   ```
   sbatch -A ACCOUNT -p PARTITION src/batches/batches.sbatch --staging-dir S --ids IDS.csv --input-dir INPUT_DIR --disk-budget 20T
   ```

   The job runs batches until no sample is `eligible`. With `--n-batches N`, it stops after N batches.
4. Monitor the samples:

   ```
   python3 src/batches/batches.py status --staging-dir S --ids IDS.csv --input-dir INPUT_DIR --disk-budget 20T
   python3 src/batches/batches.py sample --staging-dir S SRR123456
   ```

5. Promote the accepted samples. With `--dry-run`, `promote` only shows the counts:

   ```
   python3 src/batches/batches.py promote --staging-dir S --valid-dir gutmap/valid --logs-dir gutmap/valid-logs --dry-run
   python3 src/batches/batches.py promote --staging-dir S --valid-dir gutmap/valid --logs-dir gutmap/valid-logs
   ```

`promote` takes the samples with gate `yes` and state `done` or `failed`. If a target dir exists, no dir of that sample moves.

## Stop and continue

- To stop the job, send `scancel JOBID` (`make stop`). Slurm stops `batches.py` and Nextflow, and Nextflow cancels its tasks. The time limit and a node failure stop the job in the same way. The samples of the batch stay `running`.
- The next `run` finishes that batch first. The samples that did not end get `stopped`, and a later batch runs them again. `run --n-batches 0` only finishes the batch.
- Each batch is a new Nextflow run without `-resume`. After each batch, `run` removes the work dir of that run with `nextflow clean -f <session id>`. The session id comes from `<batch>/nextflow.log`.
- Only one `run` or `promote` at a time can use a staging dir. The `singleton` dependency of `batches.sbatch` stops a second job until the first ends. `promote` stops if a sample is `running`.

## Requirements

- Python 3.7 or later, standard library only.
- Nextflow 25.10.2 (`module load Nextflow/25.10.2` in `batches.sbatch`).

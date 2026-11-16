# gutMAP

gutMAP processes public human gut metagenomes from SRA and ENA on the Helios cluster of Cyfronet. The dataset has about 100,000 runs.

For each sample, the workflow gives a taxonomic profile, an assembly, an annotated gene catalog with gene abundances, and high-quality genomes from binning.

## Workflow

```
IDs (CSV)
  │  download.py runinfo, plan, run          kingfisher: aws-http or ena-ftp
  ▼
<input-dir>.download/downloaded/<id>/       file set checked
  │  download.py validate                    vdb-validate, gzip -t
  ▼
<input-dir>/<id>/                            input of the pipeline
  │  batches.py run                          one Nextflow run of mtg_pipeline.nf for each batch
  ▼
<staging-dir>/batch_NNNN/published/         results and logs of each batch
  │  batches.py promote                      accepted samples only
  ▼
<valid-dir>/, <logs-dir>/                    final dataset
```

- **Download** (`src/download/`): kingfisher downloads each run into a staging dir. Only samples that pass the file-set check and the content check go into the input dir. A sample gets at most 3 attempts.
- **Batches** (`src/batches/`): a batch takes samples until their disk reserves fill the disk budget. The reserve of a sample comes from its size in the NCBI runinfo.
- Each batch is a new Nextflow run without `-resume`, so a step can delete the reads of the step before it. `<staging-dir>/samples.tsv` records the state of each sample.
- If a batch stops, the next run finishes it. The samples that did not end go into a later batch.
- **Promote** is a manual step. It moves the results of the accepted samples to the valid dir, and their logs to the logs dir.
- **Job history** (`utils/collect_hjh.py`): it collects the Slurm jobs and their billing from `hpc-jobs-history`.
- Downloads and batches do not run at the same time. The `Makefile` runs all commands with the settings of `config.mk`. `make help` lists the targets.
- The READMEs in `src/download/` and `src/batches/` give the `make` steps in order.

## Pipeline

`src/pipeline/mtg_pipeline.nf` runs these steps for each sample. Three branches start after the shared steps.

| Part | Process | Tool | Result in the valid dir |
|---|---|---|---|
| Shared | `CHECK_AND_CONVERT_READS` | fasterq-dump, for `.sra` inputs only | – |
| | `FASTP_QC` | fastp | QC report, in the logs dir |
| | `REMOVE_HUMAN_READS` | bowtie2 against hg38 | – |
| | `MEGAHIT_ASSEMBLY` | MEGAHIT, contigs of 1,500 bp or more | `04_shared_assembly/` |
| Branch 1: taxonomy | `SYLPH_TAXONOMY` | sylph with GTDB r226, on the non-human reads | `branch_1/01_taxonomy/` |
| Branch 2: genes | `BAKTA_ANNOTATION` | Bakta 1.12.2 | `branch_2/01_annotation/` |
| | `CREATE_CATALOG_AND_INDEX` | cd-hit-est at 100 % identity, bwa-mem2 index | `branch_2/02_gene_catalog/` |
| | `ALIGN_AND_QUANTIFY_READS` | bwa-mem2, samtools idxstats | `branch_2/03_gene_quantification/` |
| | `CALCULATE_TPM_AND_ANNOTATE` | `bin/tpm_annotate.py` | `branch_2/04_tpm_and_annotate/` |
| Branch 3: genomes | `MAP_FOR_BINNING` | bowtie2 against the contigs, samtools | – |
| | `METABAT2_BINNING` | MetaBAT2 | – |
| | `CHECKM_QA` | CheckM2 | `branch_3/03_checkm_qa/` |
| | `FILTER_AND_ANNOTATE` | `bin/filter_bins.py`, GTDB-Tk `classify_wf` (bac120) | `branch_3/04_filter_and_annotate/` |

- **Acceptance:** a sample is accepted when its sylph profile is larger than 218 B (more than the header) and its contigs file is not empty.
- **Empty input:** a sample with no non-human reads, no contigs or no bins ends that branch without a failure.
- **High-quality bins:** completeness 90 % or more and contamination 5 % or less, from CheckM2.
- **Retry:** a failed task does not run again. The exception is `FILTER_AND_ANNOTATE`, which runs once more with 35 cores and 140 GB.
- **Slurm:** `nextflow.config` sets the account and the partitions `plgrid` and `plgrid-bigmem`. Each task is one Slurm job.

## Requirements

- Helios, with a Slurm account that can use the partitions `plgrid` and `plgrid-bigmem`.
- Nextflow 25.10.2, from `module load Nextflow/25.10.2`.
- Conda. Nextflow makes the environment of each step from bioconda, in the cache dir of `nextflow.config`.
- A conda env with kingfisher 0.5.0 and sra-tools 3.4.1 for the downloads. `config.mk` gives its path as `KF_ENV`.
- Python 3.7 or later, standard library only, for `download.py`, `batches.py` and `collect_hjh.py`.
- The databases, with their paths in `nextflow.config`:

  | Database | Used by |
  |---|---|
  | hg38 bowtie2 index | `REMOVE_HUMAN_READS` |
  | sylph GTDB r226 (`gtdb-r226-c200-dbv1.syldb`) | `SYLPH_TAXONOMY` |
  | Bakta DB v6 | `BAKTA_ANNOTATION` |
  | CheckM2 DB (`uniref100.KO.1.dmnd`) | `CHECKM_QA` |
  | GTDB-Tk release 226 | `FILTER_AND_ANNOTATE` |

- The input dir, the staging dir, the valid dir and the logs dir must be on one filesystem, because `promote` moves the dirs with `os.rename`.
- `NCBI_API_KEY` in the environment is optional. With a key, NCBI allows more requests for the runinfo.

## Outputs

**Final dataset.** Each sample has one dir in each step that gave a result:

```
<valid-dir>/04_shared_assembly/<id>/              final.contigs.fa
<valid-dir>/branch_1/01_taxonomy/<id>/            <id>.sylph_profile.tsv
<valid-dir>/branch_2/01_annotation/<id>/          <id>.faa, .ffn, .tsv, .gff3, .txt
<valid-dir>/branch_2/02_gene_catalog/<id>/        <id>_clustered_catalog.fna
<valid-dir>/branch_2/03_gene_quantification/<id>/ <id>.idxstats.txt
<valid-dir>/branch_2/04_tpm_and_annotate/<id>/    <id>_gene_quantification.tsv
<valid-dir>/branch_3/03_checkm_qa/<id>/           quality_report.tsv, DIAMOND_RESULTS.tsv
<valid-dir>/branch_3/04_filter_and_annotate/<id>/ hq_bins/, gtdbtk.bac120.summary.tsv or NO_HQ_BINS_FOUND.tsv
<logs-dir>/<step>/<id>/                           the logs of each step; the fastp reports in 02_fastp_qc
```

**Records.**

| File | Content |
|---|---|
| `<input-dir>.download/runinfo.csv` | the NCBI runinfo row of each sample |
| `<input-dir>.download/attempts.csv` | one row for each download attempt |
| `<input-dir>.download/validations.csv` | one row for each content check |
| `<staging-dir>/samples.tsv` | the state of each sample, the acceptance result and the result of each step |
| `<staging-dir>/traces.tsv` | the Nextflow traces of all batches |
| `<staging-dir>/not_planned.tsv` | the samples that are too large for the disk budget or have no size |
| `<staging-dir>/promote.csv` | one row for each dir that `promote` moved |
| `<staging-dir>/batch_NNNN/failed_logs/` | the logs of the failed tasks |
| `HJH_CSV` in `config.mk` | the Slurm jobs and their billing, from `utils/collect_hjh.py` |

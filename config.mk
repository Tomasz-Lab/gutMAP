# Settings of the Makefile. Set the values before the first run.
# Relative paths start at the launch dir. NCBI_API_KEY and NCBI_EMAIL come only from the environment.

IDS            = data/samples.csv
INPUT_DIR      = /path/to/input_mg/input_dataset
STAGING_DIR    = /path/to/output_mg/staging
VALID_DIR      = /path/to/output_mg/gutmap/valid
VALID_LOGS_DIR = /path/to/output_mg/gutmap/valid-logs
DISK_BUDGET    = 20T
KF_ENV         = /path/to/conda/env-with-kingfisher

# Slurm account and partition. Leave them empty to use SBATCH_ACCOUNT and SBATCH_PARTITION.
ACCOUNT        =
PARTITION      =

JOB_LOG_DIR    = logs
METHODS        =
HJH_CSV        = hpc-jobs-history.csv
HJH_START      =

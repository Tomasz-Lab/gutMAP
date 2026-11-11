#!/bin/bash
#SBATCH --job-name=gutmap-validate
#SBATCH --dependency=singleton
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#
# download.py validate, with the arguments of this script. KF_ENV is a conda env with kingfisher and sra-tools.
# Submit from the dir of download.py: Slurm runs a copy of this file, so $0 does not give that dir.
# A --dependency on the command line replaces the one in this file, so keep singleton in it:
#
#   KF_ENV=/path/to/env sbatch -A ACCOUNT -p PARTITION --dependency=singleton,afterany:DOWNLOAD_JOB \
#       validate.sbatch --input-dir INPUT_DIR
set -euo pipefail
: "${KF_ENV:?KF_ENV must be the conda env with kingfisher and sra-tools}"
export PATH="$KF_ENV/bin:$PATH"
exec python3 "$SLURM_SUBMIT_DIR/download.py" validate --threads "$SLURM_CPUS_PER_TASK" "$@"

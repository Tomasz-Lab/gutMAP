#!/bin/bash
#SBATCH --job-name=gutmap-download
#SBATCH --dependency=singleton
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#
# download.py run, with the arguments of this script. KF_ENV is a conda env with kingfisher and sra-tools.
# Submit from the dir of download.py: Slurm runs a copy of this file, so $0 does not give that dir.
#
#   KF_ENV=/path/to/env sbatch -A ACCOUNT -p PARTITION download.sbatch --ids samples.csv --input-dir INPUT_DIR
set -euo pipefail
: "${KF_ENV:?KF_ENV must be the conda env with kingfisher and sra-tools}"
export PATH="$KF_ENV/bin:$PATH"
exec python3 "$SLURM_SUBMIT_DIR/download.py" run "$@"

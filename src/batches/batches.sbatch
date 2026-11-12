#!/bin/bash -l
#SBATCH --job-name=gutmap-batches
#SBATCH --dependency=singleton
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=72:00:00
#SBATCH --signal=B:TERM@900
#
# batches.py run, with the arguments of this script. Nextflow runs in this job and submits the tasks.
# Submit from the launch dir of the pipeline, the dir that holds src/. Slurm runs a copy of this file,
# so $0 does not give the dir of batches.py.
#
#   sbatch -A ACCOUNT -p PARTITION src/batches/batches.sbatch --staging-dir S --ids IDS.csv \
#       --input-dir INPUT_DIR --disk-budget 20T
#
# To stop: scancel --signal=TERM --batch JOBID. Only batches.py gets the signal. It stops nextflow,
# finishes the batch and ends. A plain scancel stops all processes at the same time, and the next
# run finishes the batch.
set -eo pipefail
module load Nextflow/25.10.2
set -u
exec python3 "$SLURM_SUBMIT_DIR/src/batches/batches.py" run "$@"

#!/usr/bin/env bash
# run_batches.sh — runs mtg_pipeline.nf in batches over samplesheet
#
# Usage (run from project root):
#   bash src/run_batches.sh [--start N] [--end N] [--batch-size N]
#
# --start      First sample index to process, 1-based (default: 1)
# --end        Last sample index to process (default: all)
# --batch-size Samples per Nextflow run (default: 500)
#
# Example: process samples 1001-5000 in batches of 500
#   bash src/run_batches.sh --start 1001 --end 5000 --batch-size 500

set -euo pipefail

PROJ="/net/scratch/hscra/plgrid/plgpkica/metagenome_proj"
FULL_CSV="$PROJ/data/helios_greater_than_500MB_sorted.csv"
PIPELINE="$PROJ/src/mtg_pipeline_batch_bakta.nf"
WORK_DIR="/net/storage/pr3/plgrid/plgggutmap/nf_work"
BATCH_DIR="$PROJ/data/greater_than_500MB_sorted"
LOG_DIR="$PROJ/results/logs/greater_than_500MB_sorted"
OUTPUT_BASE="/net/storage/pr3/plgrid/plgggutmap/output_mg/greater_than_500MB_sorted"
INPUT_BASE="/net/storage/pr3/plgrid/plgggutmap/input_mg/input_1-15000"

TOTAL_SAMPLES=$(( $(wc -l < "$FULL_CSV") - 1 ))  # subtract header row

# --- Defaults ---
BATCH_SIZE=500
START=7501
END=9000

# --- Parse args ---
while [[ $# -gt 0 ]]; do
    case $1 in
        --start)      START="$2";      shift 2 ;;
        --end)        END="$2";        shift 2 ;;
        --batch-size) BATCH_SIZE="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# --- Validate ---
if [ "$START" -lt 1 ] || [ "$END" -gt "$TOTAL_SAMPLES" ] || [ "$START" -gt "$END" ]; then
    echo "ERROR: invalid range --start $START --end $END (CSV has $TOTAL_SAMPLES samples)"
    exit 1
fi

mkdir -p "$BATCH_DIR" "$LOG_DIR"

echo "========================================================"
echo "  Batch runner started : $(date)"
echo "  Full CSV             : $FULL_CSV"
echo "  Total samples in CSV : $TOTAL_SAMPLES"
echo "  Processing range     : $START - $END"
echo "  Batch size           : $BATCH_SIZE"
echo "  Output base dir      : $OUTPUT_BASE"
echo "  Work dir (cleaned)   : $WORK_DIR"
echo "========================================================"

BATCH_NUM=1
CURRENT=$START

# Run Nextflow
module load Nextflow/25.10.2 2>/dev/null || true

while [ "$CURRENT" -le "$END" ]; do
    BATCH_END=$(( CURRENT + BATCH_SIZE - 1 ))
    [ "$BATCH_END" -gt "$END" ] && BATCH_END="$END"

    BATCH_NAME="batch_s${CURRENT}-${BATCH_END}"
    BATCH_CSV="$BATCH_DIR/${BATCH_NAME}.csv"
    NF_LOG="$LOG_DIR/${BATCH_NAME}.log"
    BATCH_OUTDIR="$OUTPUT_BASE/$BATCH_NAME"

    # Build batch CSV: header + slice
    # Sample index N in the data = line N+1 in the file (line 1 is "sra_id" header)
    echo "sra_id" > "$BATCH_CSV"
    sed -n "$(( CURRENT + 1 )),$(( BATCH_END + 1 ))p" "$FULL_CSV" >> "$BATCH_CSV"

    ACTUAL=$(( $(wc -l < "$BATCH_CSV") - 1 ))

    echo ""
    echo "--- Batch $BATCH_NUM | samples $CURRENT-$BATCH_END ($ACTUAL IDs) | $(date '+%Y-%m-%d %H:%M:%S') ---"
    echo "    Batch CSV  : $BATCH_CSV"
    echo "    Output dir : $BATCH_OUTDIR"
    echo "    NF log     : $NF_LOG"

    nextflow \
        -log "$NF_LOG" \
        run "$PIPELINE" \
        --input "$BATCH_CSV" \
        --indir "$INPUT_BASE" \
        --outdir "$BATCH_OUTDIR" \
        --limit "$BATCH_SIZE"
    NF_EXIT=$?

    if [ $NF_EXIT -ne 0 ]; then
        echo "[WARN] Nextflow exited with code $NF_EXIT for batch $BATCH_NUM — continuing to next batch"
    else
        echo "[OK] Batch $BATCH_NUM finished successfully"
    fi

    # Clean up work dir to free scratch space
    if [ -d "$WORK_DIR" ]; then
        echo "Cleaning up $WORK_DIR ..."
        rm -rf "$WORK_DIR"
        nextflow clean -f || true
        echo "Done."
    fi

    CURRENT=$(( BATCH_END + 1 ))
    BATCH_NUM=$(( BATCH_NUM + 1 ))
done

echo ""
echo "========================================================"
echo "  All batches complete : $(date)"
echo "  Batches run          : $(( BATCH_NUM - 1 ))"
echo "========================================================"

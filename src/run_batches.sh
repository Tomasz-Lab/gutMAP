#!/usr/bin/env bash
# run_batches.sh — run N Nextflow batches back-to-back, each over the next
# smallest unprocessed sra_ids (auto mode); size workdirs after every batch.
#
# Usage:
#   bash src/run_batches.sh [--batch-size N] [--n-batches N] [--csv path.csv]
#
# --batch-size N   Ids per batch when auto-selecting (default: 500).
# --n-batches  N   Number of consecutive batches to run (default: 1).
#                  Exits early if the master CSV is exhausted before N is hit.
# --csv path.csv   Skip auto-selection; use this file verbatim. File must have
#                  a `sra_id` header. The CSV is copied into $BATCH_DIR for the
#                  audit trail so future auto runs will skip its ids. Mutually
#                  exclusive with --n-batches > 1.
#
# Auto-select picks the next N smallest unprocessed ids from $FULL_CSV (sorted
# ascending by Bytes). "Processed" = present in any of:
#   - $GUTMAP_BASE/{valid,retry,invalid}/02_fastp_qc/<sra_id>/   (classified)
#   - $OUTPUT_BASE/*/published/02_fastp_qc/<sra_id>/             (pending classification)
#   - $BATCH_DIR/*.csv                                           (already submitted)
#
# Failed batches are NOT auto-retried; use --csv for targeted retries.

set -euo pipefail

PROJ="/net/scratch/hscra/plgrid/plgpkica/metagenome_proj"
FULL_CSV="$PROJ/data/input_dataset.csv"
PIPELINE="$PROJ/src/mtg_pipeline.nf"
WORK_DIR="/net/storage/pr3/plgrid/plgggutmap/nf_work"
BATCH_DIR="$PROJ/data/input_dataset"
LOG_DIR="$PROJ/results/logs/input_dataset"
INPUT_BASE="/net/storage/pr3/plgrid/plgggutmap/input_mg/input_dataset"
OUTPUT_BASE="/net/storage/pr3/plgrid/plgggutmap/output_mg/output_dataset"
GUTMAP_BASE="/net/storage/pr3/plgrid/plgggutmap/output_mg/gutmap"
CANONICAL_STEP="02_fastp_qc"   # first per-sra published dir

# --- Defaults ---
BATCH_SIZE=750
N_BATCHES=7
CSV_OVERRIDE=""

# --- Parse args ---
while [[ $# -gt 0 ]]; do
    case $1 in
        --batch-size) BATCH_SIZE="$2"; shift 2 ;;
        --n-batches)  N_BATCHES="$2";  shift 2 ;;
        --csv)        CSV_OVERRIDE="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [ -n "$CSV_OVERRIDE" ] && [ "$N_BATCHES" -gt 1 ]; then
    echo "ERROR: --csv is mutually exclusive with --n-batches > 1" >&2
    exit 1
fi

mkdir -p "$BATCH_DIR" "$LOG_DIR" "$PROJ/data"

# Reusable temp files (overwritten each iteration). PID suffix avoids
# collisions if two run_batches.sh ever overlap on the same node.
SKIP="$PROJ/data/.skip_set.$$"
DOWNLOADED="$PROJ/data/.downloaded.$$"
trap 'rm -f "$SKIP" "$DOWNLOADED"' EXIT

module load Nextflow/25.10.2 2>/dev/null || true

# ============================================================
# Main loop: N consecutive batches (or until exhausted)
# ============================================================
for (( i=1; i<=N_BATCHES; i++ )); do
    echo
    echo "######################  Batch $i / $N_BATCHES  ######################"

    # --- 1. Pick the batch CSV (override path OR auto-select) ---
    if [ -n "$CSV_OVERRIDE" ]; then
        [ -f "$CSV_OVERRIDE" ] || { echo "ERROR: not a file: $CSV_OVERRIDE"; exit 1; }
        BATCH_NAME="batch_override_$(date +%Y%m%d_%H%M%S)"
        BATCH_CSV="$BATCH_DIR/${BATCH_NAME}.csv"
        cp "$CSV_OVERRIDE" "$BATCH_CSV"
        ACTUAL=$(( $(wc -l < "$BATCH_CSV") - 1 ))
        echo "[override] $ACTUAL ids -> $BATCH_CSV"
    else
        # Eligibility filter: what's currently downloaded on disk. New downloads
        # added between runs flow in automatically here.
        if [ -d "$INPUT_BASE" ]; then
            ls -1 "$INPUT_BASE" | sort -u > "$DOWNLOADED"
        else
            : > "$DOWNLOADED"
        fi

        {
            # Classified done: gutmap/{valid,invalid}/02_fastp_qc/<sra_id>/
            for status in valid invalid; do
                d="$GUTMAP_BASE/$status/$CANONICAL_STEP"
                [ -d "$d" ] && ls -1 "$d"
            done
            # Pending classification: output_dataset/<batch>/published/02_fastp_qc/<sra_id>/
            # (mid-flight or awaiting manual mv to gutmap/).
            for d in "$OUTPUT_BASE"/*/published/"$CANONICAL_STEP"; do
                [ -d "$d" ] && ls -1 "$d"
            done
            :
        } | sort -u > "$SKIP"

        echo "[auto] downloaded: $(wc -l < "$DOWNLOADED")  skip set: $(wc -l < "$SKIP")"

        # Next sequence number (batch_NNNN.csv only; override CSVs ignored)
        SEQ=$(find "$BATCH_DIR" -maxdepth 1 -name 'batch_[0-9]*.csv' | wc -l)
        BATCH_NAME=$(printf 'batch_%04d' "$((SEQ + 1))")
        BATCH_CSV="$BATCH_DIR/${BATCH_NAME}.csv"

        # input_dataset.csv is sorted ascending by Bytes -> walking top-to-bottom and
        # filtering by (downloaded AND not-in-skip) yields smallest-unprocessed-first
        # among samples whose input is actually on disk.
        # Three input files; FILENAME idiom disambiguates which set each row builds.
        # Single awk pass exits on LIMIT (no `head -n` -> no SIGPIPE).
        echo "sra_id" > "$BATCH_CSV"
        awk -F, -v LIMIT="$BATCH_SIZE" -v skipf="$SKIP" -v dlf="$DOWNLOADED" '
            FILENAME == skipf { skip[$0] = 1; next }
            FILENAME == dlf   { dl[$0]   = 1; next }
            FNR == 1          { next }              # drop $FULL_CSV header
            {
                id = $1
                if ((id in dl) && !(id in skip)) {
                    print id
                    if (++n >= LIMIT) exit
                }
            }
        ' "$SKIP" "$DOWNLOADED" "$FULL_CSV" >> "$BATCH_CSV"

        ACTUAL=$(( $(wc -l < "$BATCH_CSV") - 1 ))
        if [ "$ACTUAL" -eq 0 ]; then
            echo "[auto] no unprocessed sra_ids left — stopping after $((i - 1)) batch(es)"
            rm -f "$BATCH_CSV"
            break
        fi
        echo "[auto] $BATCH_NAME: $ACTUAL ids -> $BATCH_CSV"
    fi

    NF_LOG="$LOG_DIR/${BATCH_NAME}.log"
    BATCH_OUTDIR="$OUTPUT_BASE/$BATCH_NAME"

    echo "========================================================"
    echo "  Batch        : $BATCH_NAME"
    echo "  IDs          : $ACTUAL"
    echo "  Started      : $(date '+%Y-%m-%d %H:%M:%S')"
    echo "  Batch CSV    : $BATCH_CSV"
    echo "  NF log       : $NF_LOG"
    echo "  Output dir   : $BATCH_OUTDIR"
    echo "========================================================"

    # --- 2. Run Nextflow ---
    set +e
    nextflow \
        -log "$NF_LOG" \
        run "$PIPELINE" \
        --input "$BATCH_CSV" \
        --indir "$INPUT_BASE" \
        --outdir "$BATCH_OUTDIR" \
        --limit "$BATCH_SIZE"
    NF_EXIT=$?
    set -e

    if [ $NF_EXIT -ne 0 ]; then
        echo "[WARN] Nextflow exited with code $NF_EXIT for $BATCH_NAME — continuing to sizing"
    else
        echo "[OK] $BATCH_NAME finished successfully"
    fi

    # --- 3. Per-step workdir sizing ---
    echo "Submitting sizing job for $BATCH_NAME (sbatch --wait) ..."
    if sbatch --wait "$PROJ/src/measure_workdirs.sbatch" "$NF_LOG" "$BATCH_NAME"; then
        echo "[OK] sizing complete for $BATCH_NAME"
    else
        echo "[WARN] sizing job failed (rc=$?) — continuing"
    fi

     # Clean up work dir to free scratch space (runs between batches if enabled)
     if [ -d "$WORK_DIR" ]; then
         echo "Cleaning up $WORK_DIR ..."
         rm -rf "$WORK_DIR"
         nextflow clean -f || true
         echo "[OK] cleanup done"
     fi
done

echo
echo "All done."

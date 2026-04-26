#!/bin/bash
#SBATCH --job-name=s2m_arr
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --array=1-61
#SBATCH --output=/home/yxu220001/work/logs/s2m_arr_%A_%a.out
#SBATCH --error=/home/yxu220001/work/logs/s2m_arr_%A_%a.err

set -euo pipefail
export PATH=/home/yxu220001/pipeline_env/bin:$PATH

INPUT_DIR=/home/yxu220001/work/pipeline_output/bowtie2_eocrc
OUT_DIR=/home/yxu220001/work/strainphlan/consensus_markers
DB=/home/yxu220001/work/databases/metaphlan4/mpa_vJan25_CHOCOPhlAnSGB_202503.pkl

mkdir -p "$OUT_DIR"

# Get the i-th SAM file (sorted alphabetically for deterministic mapping)
SAM_FILE=$(ls $INPUT_DIR/*.sam.bz2 | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
SAMPLE=$(basename "$SAM_FILE" .sam.bz2)

echo "=== Task $SLURM_ARRAY_TASK_ID: $SAMPLE ==="
echo "Start: $(date)"
echo "Input: $SAM_FILE"

sample2markers.py \
    -i "$SAM_FILE" \
    -o "$OUT_DIR" \
    -f bz2 \
    --database "$DB" \
    --clades t__SGB6653 t__SGB6649 t__SGB7295 t__SGB748 t__SGB5842 t__SGB6007 t__SGB6031 \
    -n 8

echo "=== $SAMPLE COMPLETE: $(date) ==="
ls -la "$OUT_DIR/${SAMPLE}"* 2>/dev/null || echo "(no output for this sample - likely no target species detected)"

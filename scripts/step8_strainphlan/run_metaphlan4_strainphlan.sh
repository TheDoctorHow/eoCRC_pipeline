#!/bin/bash
#SBATCH --job-name=mp4_eocrc
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --array=1-61
#SBATCH --output=/home/yxu220001/work/logs/mp4_eocrc_%A_%a.out
#SBATCH --error=/home/yxu220001/work/logs/mp4_eocrc_%A_%a.err

set -euo pipefail

# Activate environment
export PATH=/home/yxu220001/pipeline_env/bin:$PATH

# Parameters
MANIFEST=/home/yxu220001/work/strainphlan/samples_paired.tsv
METAPHLAN_DB=/home/yxu220001/work/databases/metaphlan4
OUTDIR=/home/yxu220001/work/pipeline_output/bowtie2_eocrc
THREADS=16

# Row (skip header): array task ID i means line i+1
ROW=$(awk -v i="$SLURM_ARRAY_TASK_ID" 'NR==i+1' "$MANIFEST")
SAMPLE=$(echo "$ROW" | cut -f1)
STUDY=$(echo "$ROW" | cut -f2)
N_RUNS=$(echo "$ROW" | cut -f3)
R1_LIST=$(echo "$ROW" | cut -f4)
R2_LIST=$(echo "$ROW" | cut -f5)

echo "=== Task $SLURM_ARRAY_TASK_ID: $SAMPLE ($STUDY, $N_RUNS runs) ==="
echo "Start: $(date)"

# Scratch setup
SCRATCH=/scratch/juno/yxu220001/mp4_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCRATCH"
trap "rm -rf $SCRATCH" EXIT

# Concatenate runs
R1_FILES=$(echo "$R1_LIST" | tr ';' ' ')
R2_FILES=$(echo "$R2_LIST" | tr ';' ' ')

R1_CAT="$SCRATCH/${SAMPLE}_R1.fastq.gz"
R2_CAT="$SCRATCH/${SAMPLE}_R2.fastq.gz"

echo "--- Concatenating all reads (R1+R2 as single-end for StrainPhlAn) ---"
zcat $R1_FILES $R2_FILES | pigz -p $THREADS -c > "$R1_CAT"

echo "Concat sizes:"
ls -lh "$R1_CAT"

# Run MetaPhlAn4 with bowtie2out for StrainPhlAn
echo "--- MetaPhlAn4 start: $(date) ---"
metaphlan \
    "$R1_CAT" \
    --input_type fastq \
    --db_dir "$METAPHLAN_DB" \
    --nproc $THREADS \
    --force \
    --mapout "$SCRATCH/${SAMPLE}.mapout.bz2" \
    -s "$SCRATCH/${SAMPLE}.sam.bz2" \
    --output_file "$SCRATCH/${SAMPLE}_profile.tsv"
echo "--- MetaPhlAn4 done: $(date) ---"

# Copy outputs to work
cp "$SCRATCH/${SAMPLE}.sam.bz2" "$OUTDIR/"
cp "$SCRATCH/${SAMPLE}_profile.tsv" "$OUTDIR/"

echo "=== $SAMPLE COMPLETE: $(date) ==="
ls -lh "$OUTDIR/${SAMPLE}"*

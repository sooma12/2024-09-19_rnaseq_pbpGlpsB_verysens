#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=alignRNA_bowtie
#SBATCH --time=08:00:00
#SBATCH --array=<>
#SBATCH --ntasks=<>
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8
#SBATCH --output=<>
#SBATCH --error=<>
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<>

## Usage: sbatch 3_sbatch_array_bowtie2_align.sh
echo "Loading tools"
module load bowtie/2.5.2

# Load config file
source ./config.cfg

echo "Bowtie2 reference genome index basename: $BT2_OUT_BASE"
echo "sample sheet located at $SAMPLE_SHEET_PATH"
echo "alignment output directory containing .bam files: ${MAPPED_DIR}"

mkdir -p $MAPPED_DIR

# Array job!  Used my sample sheet technique from 2023-11-13 breseq for this.
# NOTE!!! sample sheet prep is moved to separate script.
# Run that script before this one!

# sed and awk read through the sample sheet and grab each whitespace-separated value
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p $SAMPLE_SHEET_PATH |  awk '{print $1}')
r1=$(sed -n "$SLURM_ARRAY_TASK_ID"p $SAMPLE_SHEET_PATH |  awk '{print $2}')
r2=$(sed -n "$SLURM_ARRAY_TASK_ID"p $SAMPLE_SHEET_PATH |  awk '{print $3}')

echo "Running Bowtie2 on files $r1 and $r2"

# Bowtie2 in paired end mode
bowtie2 --local -p 8 -x $BT2_OUT_BASE -q -1 $r1 -2 $r2 -S $MAPPED_DIR/$name.sam

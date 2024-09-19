#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=samtools_to_sorted_bam
#SBATCH --time=04:00:00
#SBATCH --array=<>
#SBATCH --ntasks=<>
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --output=<>
#SBATCH --error=<>
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<>

echo "Loading tools"
module load samtools/1.19.2

source ./config.cfg

echo "Looking for .sam files and outputting sorted bams in $MAPPED_DIR"

# Get array of sam files
# shellcheck disable=SC2207
sams_array=($(ls -d ${MAPPED_DIR}/*.sam))

# Get specific file for this array task
current_file=${sams_array[$SLURM_ARRAY_TASK_ID-1]}

current_name=$(basename "$current_file")
current_name_no_ext="${current_name%.*}"

echo "converting .sam to .bam"
samtools view -bS "${current_file}" > ${MAPPED_DIR}/${current_name_no_ext}.bam
echo "sorting bam file"
samtools sort ${MAPPED_DIR}/${current_name_no_ext}.bam -o "${MAPPED_DIR}/${current_name_no_ext}"_sorted.bam
echo "analyzing alignment statistics"
samtools stats "${MAPPED_DIR}/${current_name_no_ext}"_sorted.bam > "${MAPPED_DIR}/${current_name_no_ext}"_sorted.bam.stats

echo "cleaning up: moving .sam and unsorted .bam files to $MAPPED_DIR/intermediate_files"
mkdir -p $MAPPED_DIR/intermediate_files
echo "moving ${current_file}"
mv ${current_file} $MAPPED_DIR/intermediate_files
echo "moving ${current_file}${current_name_no_ext}.bam"
mv $MAPPED_DIR/"${current_name_no_ext}".bam $MAPPED_DIR/intermediate_files

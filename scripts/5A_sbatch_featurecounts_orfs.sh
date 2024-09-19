#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=featureCounts
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=<>
#SBATCH --error=<>
#SBATCH --mail-type=END
#SBATCH --mail-user=<>

echo "Loading environment and tools"
module load anaconda3/2021.05
eval "$(conda shell.bash hook)"
conda activate /work/geisingerlab/conda_env/subread

source ./config.cfg
echo "featureCounts file name: $ORF_COUNTS_FILE found in $COUNTS_OUTDIR"
echo "genome GTF reference file: $GENOME_GTF"
echo ".bam input files were found in: $MAPPED_DIR"

mkdir -p $COUNTS_OUTDIR

# Run featureCounts on all BAM files from STAR
# -t flag specifies features in column 3 of GTF to count; default is exon.
featureCounts \
-a $GENOME_GTF \
-o $COUNTS_OUTDIR/$ORF_COUNTS_FILE \
-p \
--countReadPairs \
-t transcript \
$MAPPED_DIR/*.bam

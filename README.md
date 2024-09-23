# 2024-09-17_pbpGlpsB_very-sensitive

1. Link to files in `2024-01_rnaseq_pbpGlpsB`
`for file in /work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/input/fastq/
do
ln -s $file input/
done`

2. Bowtie build
Used script 1; default Bowtie-build parameters

3. Sample sheet
Script 2

4. Bowtie2 alignment

Ran bowtie2 aligner with following arguments:
--local
--very-sensitive-local
-p 8 
-x $BT2_OUT_BASE 
-q 
-1 $r1 
-2 $r2 
-S $MAPPED_DIR/$name.sam

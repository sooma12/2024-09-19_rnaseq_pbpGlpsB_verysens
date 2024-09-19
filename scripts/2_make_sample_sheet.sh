#!/bin/bash
# 2_make_sample_sheet.sh
# Makes a sample_sheet.txt containing sample ID and R1 and R2 filepaths
# Assumes each sample file is in the format: WT_1_S84_trimmed_R2.fastq
# Script uses a split on "_" to grab the sample ID, e.g. WT_1.  Must modify this if sample file names are different!
# Example output line:  WT_1 /path/to/WT_1_R1.fastq /path/to/WT_1_R2.fastq

source ./config.cfg

# Create .list files with R1 and R2 fastqs.  Sort will put them in same orders, assuming files are paired
find ${FASTQDIR_UNTRIMMED} -maxdepth 1 -name "*.fastq" | grep -e "_R1" | sort > R1.list
find ${FASTQDIR_UNTRIMMED} -maxdepth 1 -name "*.fastq" | grep -e "_R2" | sort > R2.list

if [ -f "${SAMPLE_SHEET_PATH}" ] ; then
  rm "${SAMPLE_SHEET_PATH}"
fi

# make sample sheet from R1 and R2 files.  Format on each line looks like (space separated):
# WT_1 /path/to/WT_1_R1.fastq /path/to/WT_1_R2.fastq
# from sample sheet, we can access individual items from each line with e.g. `awk '{print $3}' sample_sheet.txt`

paste R1.list R2.list | while read R1 R2 ;
do
    outdir_root=$(basename ${R2} | cut -f1,2 -d"_")
    sample_line="${outdir_root} ${R1} ${R2}"
    echo "${sample_line}" >> $SAMPLE_SHEET_PATH
done

rm R1.list R2.list

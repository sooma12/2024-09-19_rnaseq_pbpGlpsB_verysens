# analyze_pca_pbpGlpsB_2024-09-24.R

# https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html#quick-start-gene-expression-omnibus-geo

#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')

#BiocManager::install('PCAtools')

# remotes::install_version("matrixStats", version="1.1.0") # restart your session and run previous scripts

library('PCAtools')
library(dplyr)
library(tidyverse)
library('DESeq2')


######
### PCAs for ORFs
######

## Read in raw data and prepare feature count table
feature_count <- read.table("~/Documents/geisinger_lab_research/bioinformatics_in_acinetobacter/rnaSeq/2024-09-17_pbpGlpsB_very-sensitive/data/counts_orfs.txt", header=TRUE, row.names = 1)
# Grab the count data
data <- feature_count[,6:14]
# view(data)

## Set up column names in feature count data AND prep metadata table with strains and conditions
# Trim column names down to just the sample IDs
column_names <- colnames(data)
column_names <- sub("X.work.geisingerlab.Mark.rnaSeq.2024.09.19_rnaseq_pbpGlpsB_verysens.data.mapped.", "", column_names)
column_names <- sub("_sorted.bam", "", column_names)
column_names
colnames(data) <- column_names
# Use regex to get condition (mutation) from strain IDs
# ".+?(?=_)":   .+?  means "match any character (.) any number of times (+?)"
# (?=_): a positive lookahead: find where the character - is matched
conditions <- str_extract(column_names, ".+?(?=_)")
# Use column_names and conditions to make metadata table
meta <- data.frame(column_names, conditions)
colnames(meta) <- c('id', 'condition')
meta  # Verify that IDs and conditions are as expected
# Write out metadata table (save!!)
#meta_filepath <- '/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data'
#meta_file <- file.path(meta_filepath, 'metadata.txt')
#write.table(meta, meta_file, row.names=FALSE)  # TODO: Come back to this. Ok to have "id" label, or needs to be blank as in examples?  May want to change colnames before this.

## Load data, pre-filter low-count genes, and relevel to set WT as the reference
des_data <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ condition)
smallestGroupSize <- 3  # should be size of smallest group; I did 3 replicates
keep <- rowSums(counts(des_data) >= 10) >= smallestGroupSize  # keep data where count is >10 in all 3 samples
des_data <- des_data[keep,]
# relevel dds$condition to set WT as the reference
des_data$condition <- relevel(des_data$condition, ref = "WT")

dds <- DESeq(des_data)

vst <- assay(vst(dds))

pca_metadata <- data.frame(row.names = colnames(vst))
p <- pca(vst, metadata = pca_metadata, removeVar = 0.1)

screeplot(p, axisLabSize = 18, titleLabSize = 22)

# Biplot: PCA plot with PC1 and 2, plus arrows indicating how much factors influence that PC

biplot(p, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)

ggsave(filename = "PCA_biplot_ORFs_d_pbpG_lpsB_2024-09-24.tiff", width = 10, height = 10)


plotloadings(p, labSize = 3)

######
### PCAs for sRNAs
######

## Read in raw data and prepare feature count table
feature_count_srna <- read.table("~/Documents/geisinger_lab_research/bioinformatics_in_acinetobacter/rnaSeq/2024-09-17_pbpGlpsB_very-sensitive/data/counts_srnas.txt", header=TRUE, row.names = 1)
# Grab the count data
data_srna <- feature_count_srna[,6:14]
# view(data)

## Set up column names in feature count data AND prep metadata table with strains and conditions
# Trim column names down to just the sample IDs
column_names_srna <- colnames(data_srna)
column_names_srna <- sub("X.work.geisingerlab.Mark.rnaSeq.2024.09.19_rnaseq_pbpGlpsB_verysens.data.mapped.", "", column_names_srna)
column_names_srna <- sub("_sorted.bam", "", column_names_srna)
column_names_srna
colnames(data_srna) <- column_names_srna
# Use regex to get condition (mutation) from strain IDs
# ".+?(?=_)":   .+?  means "match any character (.) any number of times (+?)"
# (?=_): a positive lookahead: find where the character - is matched
conditions_srna <- str_extract(column_names_srna, ".+?(?=_)")
# Use column_names and conditions to make metadata table
meta_srna <- data.frame(column_names_srna, conditions_srna)
colnames(meta_srna) <- c('id', 'condition')
meta_srna  # Verify that IDs and conditions are as expected
# Write out metadata table (save!!)
#meta_filepath <- '/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data'
#meta_file <- file.path(meta_filepath, 'metadata.txt')
#write.table(meta, meta_file, row.names=FALSE)  # TODO: Come back to this. Ok to have "id" label, or needs to be blank as in examples?  May want to change colnames before this.

## Load data, pre-filter low-count genes, and relevel to set WT as the reference
des_data_srna <- DESeqDataSetFromMatrix(countData = data_srna, colData = meta_srna, design = ~ condition)
smallestGroupSize <- 3  # should be size of smallest group; I did 3 replicates
keep_srna <- rowSums(counts(des_data_srna) >= 10) >= smallestGroupSize  # keep data where count is >10 in all 3 samples
des_data_srna <- des_data_srna[keep_srna,]
# relevel dds$condition to set WT as the reference
des_data_srna$condition <- relevel(des_data_srna$condition, ref = "WT")

dds_srna <- DESeq(des_data_srna)

vst_srna <- assay(vst(dds_srna, nsub = sum( rowMeans( counts(dds_srna, normalized=TRUE)) > 5 )))

pca_metadata_srna <- data.frame(row.names = colnames(vst_srna))
p_srna <- pca(vst_srna, metadata = pca_metadata_srna, removeVar = 0.1)

screeplot(p_srna, axisLabSize = 18, titleLabSize = 22)

# Biplot: PCA plot with PC1 and 2, plus arrows indicating how much factors influence that PC

biplot(p_srna, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)

ggsave(filename = "PCA_biplot_srnas_d_pbpG_lpsB_2024-09-24.tiff", width = 10, height = 10)


plotloadings(p, labSize = 3)

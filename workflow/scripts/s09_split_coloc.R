#!/usr/bin/Rscript

library(tidyr)
library(data.table)


coloc_path <- snakemake@input
ofile_path <- snakemake@output


# load coloc results
coloc_new <- data.table::fread(coloc_path)

# subset to "significant" coloc (PP4>0.8)
coloc_new_subset <- coloc_new[coloc_new$PP.H4.abf > 0.8, ]

# create list of snps for LD check
ld_pairs <- data.frame(
  SNP_A = coloc_new_subset$top_cond_a,
  SNP_B = coloc_new_subset$top_cond_b
  )

dim(coloc_new)
dim(coloc_new_subset)
dim(ld_pairs)
max(nchar(ld_pairs$SNP_A))
max(nchar(ld_pairs$SNP_B))

# split in case of multiple snps
ld_pairs <- separate_rows(ld_pairs, SNP_A, sep="; ")
ld_pairs <- separate_rows(ld_pairs, SNP_B, sep="; ")

# check splitted dataset
dim(ld_pairs)
max(nchar(ld_pairs$SNP_A))
max(nchar(ld_pairs$SNP_B))

# keep only unique combinations
ld_pairs <- ld_pairs %>% distinct()

dim(ld_pairs)

# save 
data.table::fwrite(ld_pairs, file = ofile_path, row.names = FALSE, quote = FALSE, sep = "\t")

#!/usr/bin/Rscript

# Required inputs to extract conditional data
# 1. identifier of protein sequence
# 2. identifier of target SNP
# 3. genomic region coordinates chr_start_end

# For more information, please refer to the issue #26 of project page on github.

library(tidyverse)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)

# columns of the lifted file
vcf_header <- c("CHROM_38", "POS_38", "SNP")


# path to MVP results
path_mvp    <- args[1] #"/scratch/dariush.ghasemi/projects/pqtl_pipeline_finemap/results/meta_hole_3M/tidy_mvp.tsv"
path_lifted <- args[2] #"/scratch/dariush.ghasemi/projects/pqtl_pipeline_finemap/results/meta_hole_3M/VCF_lifted/seq.3314.74_10_117869699_118104520_target_10:117966436:A:G.txt"
path_base   <- args[3]


# separate data file name useful for extracting SNP from RDS
lifted_name <- basename(path_lifted)
seq_name <- lifted_name %>% str_extract("seq.(\\d)+.(\\d)+")
loc_name <- lifted_name %>% str_remove("_target.*") %>% str_remove("seq.(\\d)+.(\\d)+_")
snp_name <- lifted_name %>% str_remove(".*_target_") %>% str_remove(".txt")


# output filename
out_name <- paste0(path_base, "/conditional_data_", seq_name, "_locus_", loc_name, "_target_", str_replace_all(snp_name,":","_"), ".csv")


#----------------------------------------#
#-----       Prepare MVP file       -----
#----------------------------------------#

# read MVP coloc results
mvp <- data.table::fread(path_mvp) %>%
  # access to columns requested by MVP such as gene, uniprot, etc.
  dplyr::filter(
    seqid == seq_name,
    locus == loc_name,
    SNP == snp_name
  )


#----------------------------------------#
#-----       Extract RDS file       -----
#----------------------------------------#

# 1. read RDS from path
# 2. extract conditional data for given target SNP
# 3. append MVP-desired columns to conditional data
grep_conditional <- function(cond_file, target){
  
  # read conditional data of corresponding locus
  cond_data <- readRDS(cond_file)
  
  # Access conditional data by target SNP identifier
  conditional_target <- cond_data$results[[target]] %>%
    # retain unavailable columns requested by MVP like MAF
    dplyr::mutate(
      a1  = str_extract(SNP, "([A-Z])+"),
      a2  = str_extract(SNP, "([A-Z])+$"),
      EA  = if_else(a1 == refA, a1, a2),
      NEA = if_else(a1 != refA, a1, a2),
      EAF = if_else(a1 == refA, freq_geno, 1 - freq_geno), # align EAF
      MAF = if_else(EAF <= 0.50, EAF, 1 - EAF) # compute minor allele frequency
    ) %>%     # remove unwanted columns by MVP
    dplyr::select(
      SNP, Chr, bp, EA, NEA, MAF, b, se, n, mlog10pC
    ) %>%
    # add requested columns by MVP
    dplyr::mutate(
      seqid  = mvp$seqid,
      locus  = mvp$locus,
      TISSUE = mvp$TISSUE,
      GENE_NAME = mvp$GENE_NAME,
      UNIPROT   = mvp$UNIPROT,
      PROTEIN_NAME = mvp$PROTEIN_NAME
      )
  
  return(conditional_target)
}


# iterate function to save VCF for each input target
cond_merged <- map2_dfr(mvp$file, mvp$SNP, grep_conditional)


#----------------------------------------#
#-----     Merge RDS with Lifted    -----
#----------------------------------------#

# read lifted file 
cond_merged_b38 <- data.table::fread(
  path_lifted,
  col.names = vcf_header
  ) %>%  # merge conditional data with lifted file
  right_join(
    cond_merged,
    join_by("SNP")
  )


# save conditional data file
write.csv(cond_merged_b38, file = out_name, row.names = F, quote = F)


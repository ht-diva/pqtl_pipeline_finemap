#!/usr/bin/Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

#----------#
# taking variants file as input
sets_path <- snakemake@input
file_path <- snakemake@output[["master"]]

#--------------#
# the path where coloc SNPs are saved
exmpl_path <- as.character(sets_path[1])
pcmd <- dirname(exmpl_path)

#--------------#

coloc_pattern = "seq.(\\d+).(\\d+)/finemaping/(\\d+):(\\d+):[A-Z]:[A-Z]_locus_chr(\\d+)_(\\d+)_(\\d+)_coloc_info_table.tsv"
headers = c("phenotype_id", "credible_set", "path_rds", "path_ind_snps")

# scan COJO results for all proteins sequence
files_all <- list.files(
  path = pcmd,      # the path where we scan all files in the subdirectory
  recursive = TRUE, # to show the files in subdirectories or subfolders
  full.names = TRUE # to show full path
)

# list of info tables for the input seqids
files_coloc <- grep(coloc_pattern, files_all, value = TRUE)

#--------------#
# extract seqid from input sentinel files
input_seqid <- map_dfr(
  sets_path, function(path) {
      base_path = dirname(path)
      file_name = basename(path)
      seqid = gsub(".sentinel", "", file_name) # extract the protein sequence id and remove file format
  data.frame(base_path, file_name, seqid)
  }
)

#--------------#
# extract seqid from COJO outputs
seq_list_tbl <- tibble(files_coloc) %>% dplyr::mutate(seqid = stringr::str_extract(files_coloc, "seq.\\d+.\\d+"))

# select input seqids from all present COJO outputs
files_coloc_input <- files_coloc[seq_list_tbl$seqid %in% input_seqid$seqid]

#--------------#
# Read and merge coloc info files
combined_info <- lapply(
  files_coloc_input, 
  data.table::fread, 
  header = F, 
  sep = "\t", 
  col.names = headers
  ) %>% 
  bind_rows()

#--------------#
# extract chromosomal position from file names
combined_info <- combined_info %>% dplyr::mutate(chr = stringr::str_extract(credible_set, "(\\d+)"))

#--------------#
# save the joint results
write.table(combined_info, file = file_path, quote = F, row.names = F, sep="\t")

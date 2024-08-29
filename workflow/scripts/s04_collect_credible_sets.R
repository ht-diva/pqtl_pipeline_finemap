#!/usr/bin/Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

#----------#
# taking variants file as input
sets_path <- snakemake@input
file_path <- snakemake@output[["ofile"]]

#--------------#
# the path where COJO+cond SNPs are saved
exmpl_path <- as.character(sets_path[1])
pcmd <- dirname(exmpl_path)

#--------------#
# scan COJO results for all proteins sequence
cojo_files <- list.files(
  pattern = paste0("seq.(\\d+).(\\d+)_locus_chr(\\d+)_(\\d+)_(\\d+)_conditional_snps.tsv"),
  path = pcmd,      # the path where the independents snps file live
  recursive = FALSE, # to show the files in subdirectories or subfolders
  full.names = TRUE # to show full path
)

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
seq_list_tbl <- tibble(cojo_files) %>% mutate(seqid = str_extract(cojo_files, "seq.\\d+.\\d+"))

# select input seqids from all present COJO outputs
cojo_files_input <- cojo_files[seq_list_tbl$seqid %in% input_seqid$seqid]

#--------------#
# Merge COJO+cond SNPs characteristics
cojo_meta <- tibble(
  data.table::rbindlist(
    fill = TRUE,
    lapply(
      cojo_files_input, 
      function(x) {
        data.table::fread(x, data.table=F, fill = TRUE) %>% 
        mutate(
            study_id = stringr::str_split_fixed(basename(x), "_", 2)[,1],
            locus = stringr::str_split_fixed(basename(x), "_locus_", 2)[,2],
            locus = stringr::str_remove_all(locus, "_conditional_snps.tsv")
            )
    }
    )
  )
)

cojo_meta <- cojo_meta %>% arrange(Chr, bp)

#--------------#
# save the joint results
write.csv(cojo_meta, file = file_path, quote = F, row.names = F)

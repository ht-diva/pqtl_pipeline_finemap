#!/usr/bin/Rscript

library(data.table)
library(dplyr)

#----------#
# taking variants file as input
sets_path <- snakemake@input
file_path <- snakemake@output[["ofile"]]

#--------------#
# Load-in COJO results for all proteins sequence
credible_set_list <- lapply(
  sets_path, function(path) {
    base_path = dirname(path)
    file_name = basename(path)
    seqid = stringr::str_remove(file_name, ".sentinel") # extract the protein sequence id and remove file format
    list.files(
      pattern = paste0(seqid, "_locus_chr(\\d+)_(\\d+)_(\\d+)_ind_snps.tsv"), # pattern of output file name containing independents snps
      path = base_path, # the path where the independents snps file live
      recursive = TRUE, # to show the files in subdirectories or subfolders
      full.names = TRUE # to show full path
      )}) %>% 
  unlist()

credible_set_list

#--------------#
# Merge them
cojo_meta <- tibble(
  data.table::rbindlist(
    fill = TRUE,
    lapply(
      credible_set_list, 
      function(x) {
        data.table::fread(x, data.table=F, fill = TRUE) %>% 
        mutate(
            study_id = stringr::str_split_fixed(basename(x), "_", 2)[,1],
            locus = stringr::str_split_fixed(basename(x), "_locus_", 2)[,2],
            locus = stringr::str_remove_all(locus, "_ind_snps.tsv")
            )
    }
    )
  )
) %>% arrange(Chr)

cojo_meta

#--------------#
# save the joint results
write.csv(cojo_meta, file = file_path, quote = F, row.names = F)

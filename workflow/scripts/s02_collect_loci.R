#!/usr/bin/Rscript


library(data.table)
library(dplyr)

#----------#
# taking variants file as input
args <- commandArgs(trailingOnly = TRUE)
loci_path <- snakemake@input
file_path <- snakemake@output[["ofile"]]

#--------------#
# Merge them
loci <- tibble(
  rbindlist(
    fill = TRUE,
    lapply(
      loci_path,
      function(x) fread(x, data.table=F, fill = TRUE)
      )
    )
  ) %>% 
  arrange(chr) %>%
  filter(!is.na(chr)) %>%   #remove trait without significant signals
  filter(!(chr == 6 & !(end < 28477797 | start > 33448354)))    # remove HLA region

#--------------#
# save the joint results
write.csv(loci, file = file_path, quote = F, row.names = F)


#!/usr/bin/Rscript

library(data.table)
library(dplyr)

#----------#
# taking variants file as input
sets_path <- snakemake@input
file_path <- snakemake@output[["ofile"]]
nlrp12 <- snakemake@params[["NLRP12"]]
build <- snakemake@params[["build"]]

nlrp12 <- toupper(nlrp12)
build <- as.character(build)

#--------------#
# Load-in COJO results for all proteins sequence
credible_set_list <- lapply(
  sets_path, function(path) {
    base_path = dirname(path)
    file_name = basename(path)
    seqid = stringr::str_remove(file_name, ".sentinel") # extract the protein sequence id and remove file format
    list.files(
      pattern = paste0(seqid, "_locus_chr(\\d+)_(\\d+)_(\\d+)_conditional_snps.tsv"), # pattern of output file name containing independents snps
      path = base_path, # the path where the independents snps file live
      recursive = TRUE, # to show the files in subdirectories or subfolders
      full.names = TRUE # to show full path
      )}) %>% 
  unlist()

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
            locus = stringr::str_remove_all(locus, "_conditional_snps.tsv")
            )
    }
    )
  )
) %>%
  arrange(Chr, bp)

# take the proper boundaries respect to the build
if (nlrp12 == "YES" & build == "37") {
  #NLRP12 gene maps to 54,296,995-54,327,657 in GRCh37 coordinates.
  pos.start <- 54296995
  pos.end   <- 54327657
} else if (nlrp12 == "YES" & build == "38") {
  pos.start <- 00000000
  pos.end   <- 00000000
}

# check if the user asks to define an indicator function
if (nlrp12 == "YES") {
  cojo_meta <- cojo_meta %>%
    mutate(NLRP12 = ifelse(Chr == 19 & !(bp < pos.start | bp > pos.end), "Yes", "No"))
  cat("SNPs in the NLRP12 are indicated from the list of target SNPs.")
} else {
  cat("No filter was applied on the list of target SNPs.")
}

#--------------#
# save the joint results
write.csv(cojo_meta, file = file_path, quote = F, row.names = F)

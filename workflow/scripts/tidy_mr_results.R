#!/usr/bin/Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))


#----------------------------------------#
#-----      Inputs and Outputs      -----
#----------------------------------------#

# taking variants file as input
path_input <- snakemake@input
path_out_shrink <- snakemake@output[['unique']]
path_out_full <- snakemake@output[['annotated']]
path_base  <- snakemake@params

# ensure the class of the output paths
path_mvp <- as.character(path_input)
path_out_shrink <- as.character(path_out_shrink)
path_out_full <- as.character(path_out_full)
path_base<- as.character(path_base)


#----------------------------------------#
#-----     MVP results annotated    -----
#----------------------------------------#
# read MVP coloc results
mvp <- data.table::fread(path_mvp) %>%
  dplyr::mutate(
    seqid = str_remove_all(CHRPOS_ID, ".*_"),     # extract seqid
    locus = str_remove(locus_START_END_37, "chr"), # remove 'chr' prefix from locus
    TISSUE = "WholeBlood",           # column requested by MVP
    DATASET = "INTERVAL_CHRIS_META", # column requested by MVP
    file = paste0(path_base, "cojo/", seqid, "/conditional_data_", locus, ".rds"), # path to cojo-cond results
    filename = paste0("conditional_data_", seqid, "_locus_", locus, "_target_", str_replace_all(SNP,":","_"), ".tsv")
  )

#--------------------#
# MR results annotated by filename helping to navigate to corresponding conditional data for each target
mvp_annotated <- mvp %>%
  # remove columns added by myself for consistency with the original MR results
  dplyr::select(- seqid, - locus, - TISSUE, - DATASET, - file)
  
# save annotated results
data.table::fwrite(mvp_annotated, file = path_out_full, sep="\t")

#--------------------#
# shrink MR results to only have one signif assoc per seqid_locus_target combination
mvp_only_target <- mvp %>%
  dplyr::distinct(CHRPOS_ID, .keep_all = T) %>% # remove multiple associations for a target SNP 
  dplyr::select(
    seqid, locus, SNP, TISSUE, GENE_NAME, UniProt_ID, PROTEIN_NAME, Protein.names, DATASET, file  # columns requested by MVP
  )

# save cleansed results
data.table::fwrite(mvp_only_target, file = path_out_shrink, sep="\t")


#----------------------------------------#
#-----       VCF conversion       -----
#----------------------------------------#

# Function to write VCF file
write_vcf <- function(df, out_name) {
  
  # Create a connection to write to a file
  vcf_file <- file(out_name, "w")
  
  # Use tryCatch to ensure the file is closed properly
  tryCatch({
    # Write the VCF header
    writeLines("##fileformat=VCFv4.2", vcf_file)
    writeLines("##source=RScript", vcf_file)
    writeLines("##reference=GRCh38", vcf_file)
    writeLines("##INFO=<ID=EAF,Number=1,Type=Float,Description=Effect Allele Frequency>", vcf_file)
    writeLines("##INFO=<ID=BETA,Number=1,Type=Float,Description=Effect Size Estimate>", vcf_file)
    writeLines("##INFO=<ID=SE,Number=1,Type=Float,Description=Standard Error>", vcf_file)
    writeLines("##INFO=<ID=N,Number=1,Type=Integer,Description=Sample Size>", vcf_file)
    writeLines("##INFO=<ID=MLOG10P,Number=1,Type=Float,Description=Negative Log10 P-value>", vcf_file)
    writeLines("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", vcf_file)
    
    # Write each row as a VCF entry
    for (i in 1:nrow(df)) {
      # Extract chromosome, position, reference allele (NEA), and alternate allele (EA)
      chrom <- df$Chr[i]
      pos <- df$bp[i]
      id <- df$SNP[i]
      ref <- df$EA[i]
      alt <- df$NEA[i]
      qual <- "."
      filter <- "."
      info <- "."
         
      # Write the line to the VCF file
      line <- paste(chrom, pos, id, ref, alt, qual, filter, info, sep = "\t")
      writeLines(line, vcf_file)
    }
  },
  finally = {
    close(vcf_file) # Ensure the file connection is closed
  })
}

#--------------------#
# extract conditional results and convert them to VCF
grep_conditional <- function(cond_file, target){
  
  # read conditional data of corresponding locus
  cond_data <- readRDS(cond_file)
  
  # Access conditional data by target SNP identifier
  cond_data_target <- cond_data$results[[target]]
  
  # subset mvp result only to target SNP
  mvpt <- mvp_only_target %>% dplyr::filter(SNP == target, file == cond_file)
  
  # remove unwanted columns by MVP
  df4vcf <- cond_data_target %>%
    dplyr::mutate(
      EA  = str_extract(SNP, "([A-Z])+"), # here, effect allele is equivalent to a1
      NEA = str_extract(SNP, "([A-Z])+$") # here, non-effect allele is equivalent to a2
    ) %>%
    dplyr::select(
      SNP, Chr, bp, EA, NEA
    ) %>%
    arrange(bp) # sort SNPs for converting to VCF
  
  # name VCF file properly
  ofile_name <- paste0(path_base, "VCF/", mvpt$seqid, "_", mvpt$locus, "_target_", mvpt$SNP, ".vcf") 
  
  # Create the directory plus all necessary subdirectories
  if (!dir.exists(paste0(path_base, "VCF"))) {
    dir.create(paste0(path_base, "VCF"), recursive = TRUE)
    }

  # save VCF file
  write_vcf(df4vcf, ofile_name)
  
  # return a message showing where the VCF is saved
  return(cat("\nSave VCF here: ", ofile_name))
}


# iterate function to save VCF for each input target
map2(mvp_only_target$file, mvp_only_target$SNP, grep_conditional)

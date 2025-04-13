suppressMessages(library(optparse))

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--pipeline_path", default=NULL, help="Path where Rscript lives"),
  make_option("--coloc_guide_table", default=NULL, help="Path and filename of table listing all coloc pairwise tests to perform"),
  #make_option("--coloc_id", default=NULL, help="Id code to univocally identify the colocalisation analysis"),
  make_option("--chr_cs", default=NULL, help="Chromosomes to split the analysis in"),
  make_option("--ofile", default=NULL, help="Output filename")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Source function R functions
source(paste0(opt$pipeline_path, "funs_locus_breaker_cojo_finemap_all_at_once.R"))



###### COLOCALISTION ######

# Import table listing all coloc pairwise tests to perform
coloc_combo <- data.table::fread(opt$coloc_guide_table, data.table = F)

# Split the dataframe into a list of rows
coloc_combo_ls <- split(coloc_combo, seq(nrow(coloc_combo)))

# Perform coloc! 
coloc.full <- lapply(coloc_combo_ls, function(x){
  
# Load-in conditional dataset with precomputed lABF  
  df_cond_a <- readRDS(x$t1_path_rds)
  df_cond_b <- readRDS(x$t2_path_rds)

  # take the most significant variant and extract characteristics
  df_top_a <- df_cond_a %>% 
    slice_max(mlog10pC, n = 1) %>% # collapse values in case more than one snp exists
    summarise(
      snp  = paste(unique(snp),  collapse = "; "),
      freq = paste(unique(freq), collapse = "; "),
      freq_geno = paste(unique(freq_geno), collapse = "; "),
      mlog10pC  = paste(unique(mlog10pC),  collapse = "; ")
      )
  
  # take the most significant variant and extract characteristics
  df_top_b <- df_cond_b %>% 
    slice_max(mlog10pC, n = 1) %>% # collapse values in case more than one snp exists
    summarise(
      snp  = paste(unique(snp),  collapse = "; "),
      freq = paste(unique(freq), collapse = "; "),
      freq_geno = paste(unique(freq_geno), collapse = "; "),
      mlog10pC  = paste(unique(mlog10pC),  collapse = "; ")
    )
  
# Retrieve important info from file name
  trait_a <- paste0(x$t1_seqid)
  cojo_snp_a <- gsub(paste0(".*/", trait_a, "_(.*)_locus_.*_finemap.rds"), "\\1", x$t1_path_rds)

  trait_b <- paste0(x$t2_seqid)
  cojo_snp_b <- gsub(paste0(".*/", trait_b, "_(.*)_locus_.*_finemap.rds"), "\\1", x$t2_path_rds)
  
# Perform colocalisation for each combination of independent SNPs
#  coloc.res <- hcolo.cojo.ht(
#    df1 = df_cond_a %>% dplyr::select(snp, lABF),
#    df2 = df_cond_b %>% dplyr::select(snp, lABF)
#  )

  # perform colocalization test in standard way
  coloc.res <- coloc::coloc.abf(
    dataset1 = df_cond_a |> prepare4coloc(),
    dataset2 = df_cond_b |> prepare4coloc()
  )

  # Add top cojo SNPs and traits
  coloc.res$summary <- coloc.res$summary %>%
   t() %>% as.data.frame() %>%  # transpose and turn class into dataframe to define new features
   dplyr::mutate(
     trait_a = trait_a, 
     hit_a = cojo_snp_a, 
     top_cond_a = df_top_a$snp, 
     top_freq_a = df_top_a$freq, 
     top_freq_geno_a = df_top_a$freq_geno, 
     top_mlog10pC_a = df_top_a$mlog10pC,
     trait_b = trait_b, 
     hit_b = cojo_snp_b, 
     top_cond_b = df_top_b$snp,
     top_freq_b = df_top_b$freq,
     top_freq_geno_b = df_top_b$freq_geno,
     top_mlog10pC_b = df_top_b$mlog10pC
     )
  
  return(coloc.res)
})

# Combine coloc summary results for all tested traits pairs, add traits names, convert combined list to data frame 
only_summary_df <- lapply(coloc.full, function(x) x$summary) %>% data.table::rbindlist() %>% as.data.frame()

# first define column names, then extract seqid and independent SNP from the name of conditional data in RDS
only_summary_df <- only_summary_df %>% 
  dplyr::mutate(
    target_a = hit_a,
    target_b = hit_b,
    locus_a = hit_a,
    locus_b = hit_b,
    across(c("trait_a",   "trait_b"), ~ basename(.x)),
    across(c("target_a", "target_b"), ~ basename(.x) %>% str_remove_all("_locus_.*_finemap.rds")),
    across(c("locus_a",   "locus_b"), ~ basename(.x) %>% str_remove_all("(.*)_locus_") %>% str_remove_all("_finemap.rds"))
    ) %>%
  dplyr::select(
    trait_a, trait_b,
    locus_a, locus_b,
    target_a, target_b,
    top_cond_a, top_cond_b,
    top_freq_a, top_freq_b,
    top_freq_geno_a, top_freq_geno_b,
    top_mlog10pC_a,  top_mlog10pC_b,
    nsnps:PP.H4.abf
    ) # remove hit_a and hit_b


# save colocalization results per chromosome
data.table::fwrite(only_summary_df, file = paste0(opt$ofile), quote=F, sep="\t", na=NA)


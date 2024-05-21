suppressMessages(library(optparse))
#scipen=0

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--pipeline_path", default=NULL, help="Path where Rscript lives"),
  make_option("--chr", default=NULL, help="Locus chromosome"),
  make_option("--start", default=NULL, help="Locus starting position"),
  make_option("--end", default=NULL, help="Locus ending position"),
  make_option("--phenotype_id", default=NULL, help="Trait for which the locus boundaries have been identified - relevant in cases of molQTLs"),
  make_option("--dataset_aligned", default=NULL, help="GENOME-WIDE munged and aligned dataset file"),
  make_option("--mapping", default=NULL, help="Mapping file containing variants IDs matching with genotype bfile"),
  make_option("--p_thresh3", default=1e-04, help="Noise p-values threshold for COJO"),
  make_option("--maf", default=1e-04, help="MAF filter", metavar="character"),
  make_option("--bfile", default=NULL, help="Path and prefix name of custom LD bfiles (PLINK format .bed .bim .fam)"),
  make_option("--plink2_bin", default="/ssu/gassu/software/plink/2.00_20211217/plink2", help="Path to plink2 software"),
  make_option("--gcta_bin", default="/ssu/gassu/software/GCTA/1.94.0beta/gcta64", help="Path to GCTA software"),
  make_option("--p_thresh4", default=1e-06, help="P-value significant threshold for redefining loci boundaries post-COJO"),
  make_option("--hole", default=250000, help="Minimum pair-base distance between SNPs in different loci"),
  make_option("--cs_thresh", default=NULL, help="Percentage of credible set"),
  make_option("--study_id", default=NULL, help="Id of the study"),
  make_option("--outdir", default=NULL, help="Output directory"),
  make_option("--plink2_mem", default=NULL, help="Amount of RAM necessary for genotype extraction"),
  make_option("--plink2_threads", default=NULL, help="Number of threads for genotype extraction")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Source function R functions
source(paste0(opt$pipeline_path, "funs_locus_breaker_cojo_finemap_all_at_once.R"))

# Slightly enlarge locus by 200kb!
locus_name <- paste0(opt$chr, "_", opt$start, "_", opt$end)
opt$chr    <- as.numeric(opt$chr)
opt$start  <- as.numeric(opt$start) -100000
opt$end    <- as.numeric(opt$end) + 100000

# to avoid killing plink job, reduce resources
opt$plink2_mem <- opt$plink2_mem - 512

cat(paste("\nlocus is:", locus_name))

# GWAS input
dataset_aligned <- fread(opt$dataset_aligned, data.table=F) #%>% dplyr::filter(phenotype_id==opt$phenotype_id)

# Mapping file
pwas_map <- fread(opt$mapping, data.table=F)

cat("\nAdding original alleles from mapping to GWAS summary...")
# merge map file with GWAS results
dataset_aligned <- dataset_aligned %>%
  select(- TEST, - EXTRA) %>%
  # merge summary stats with map file. Then, SNP id matches with genotype file
  left_join(pwas_map, by = c("ID" = "PREVIOUS_ID")) %>%
  dplyr::mutate(
    snp_map = ID, # workflow needs it to report cojo results
    sdY = coloc:::sdY.est(SE, A1FREQ, N),
    type = paste0('quant')
  ) %>%
  rename(SNP = ID) #to be used by COJO to merge with genotype


cat("done.")

###############
# Perform cojo
###############

cat(paste0("\nRun COJO...\n\n"))
# Break dataframe in list of single rows
conditional.dataset <- cojo.ht(
  D=dataset_aligned,
  locus_chr=opt$chr,
  locus_start=opt$start,
  locus_end=opt$end,
  p.thresh=as.numeric(opt$p_thresh3),
  maf.thresh=as.numeric(opt$maf),
  bfile=opt$bfile,
  gcta.bin=opt$gcta_bin,
  plink.bin=opt$plink2_bin,
  plink.mem=opt$plink2_mem,
  plink.threads=opt$plink2_threads
)

#saveRDS(conditional.dataset, file=paste0("condition_data_chr15.rds"))
cat(paste0("done.\nTime to draw regional association plot..."))

# Plot conditioned GWAS sum stats
dir.create(paste0(opt$outdir), recursive = TRUE)
#pdf(paste0(opt$study_id, "_locus_chr", locus_name, "_conditioned_loci.pdf"), height=3.5*nrow(conditional.dataset$ind.snps), width=10) ### have the original loci boundaries in the name, or the slightly enlarged ones?
png(paste0(opt$study_id, "_locus_chr", locus_name, "_conditioned_loci.png"), res = 300, units = "in", height=3.5*nrow(conditional.dataset$ind.snps)+2, width=10)
plot.cojo.ht(conditional.dataset) + patchwork::plot_annotation(paste("Locus chr", locus_name))
dev.off()

cat("created!")

####################
# Locus breaker BIS
###################

cat(paste0("\nApply locus breaker and widen the locus..."))
### Repeat only on dataset that have been conditioned!!
conditional.dataset$results <- lapply(conditional.dataset$results, function(x){

  ### Check if there's any SNP at p-value lower than the set threshold. Otherwise stop here
  if(isTRUE(any(x %>% pull(pC) > -log10(opt$p_thresh4)))){
    new_bounds <- locus.breaker.p(
      x,
      p.sig    = as.numeric(opt$p_thresh4),
      p.limit  = as.numeric(opt$p_thresh3),
      hole.size= opt$hole,
      p.label  = "pC",
      chr.label= "Chr",
      pos.label= "bp")

    # Slightly enlarge locus by 200kb!
    new_bounds <- new_bounds %>% dplyr::mutate(start=as.numeric(start)-100000, end=as.numeric(end)+100000)

    # Remove SNPs not included in loci boundaries
    x %>% filter(bp >= new_bounds$start & bp <= new_bounds$end)
  }
})

cat(paste0("done."))

## Remove eventually empty dataframes (caused by p_thresh4 filter)
conditional.dataset$results <- conditional.dataset$results %>% discard(is.null)

#saveRDS(conditional.dataset, file=paste0("condition_data_chr15_up", ".rds"))



#############
# Finemapping
#############

cat("\nBegin to fine-map the locus...")
# Perform finemapping of each conditional dataset
finemap.res <- lapply(conditional.dataset$results, function(x){
  finemap.cojo(x, cs_threshold=opt$cs_thresh)
})

cat(paste0("done."))

#########################################
# Organise list of what needs to be saved
#########################################

cat("\nSaving independent signals...")
## Save independent association signals
core_file_name <- paste0(opt$study_id, "_", opt$phenotype_id)
if(opt$phenotype_id=="full") { core_file_name <- gsub("_full", "", core_file_name)}
fwrite(conditional.dataset$ind.snps, paste0(core_file_name, "_locus_chr", locus_name,"_ind_snps.tsv"), sep="\t", quote=F, na=NA)

cat("done.\nSave other lABF results...")
## Save lABF of each conditional dataset
lapply(finemap.res, function(x){
  sp_file_name <- paste0(core_file_name, "_", unique(x$cojo_snp), "_locus_chr", locus_name)
  # .rds object collecting 1) lABF, 2) beta, 3) pos for all SNPs, 3) list of SNPs in the credible set
  saveRDS(x, file=paste0(sp_file_name, "_finemap.rds")) ### cojo_snp reported in the file name   #x %>% select(-cojo_snp)
  # .tsv with 1) study id and trait (if molQTL) locus info, 2) list of SNPs in the 99% credible set, 3) path and name of correspondent .rds file and 4) path and name of correspondent ind_snps.tsv table
  #  --> append each row to a master table collecting all info from processed sum stats
  ### Idea: create guidelines for generating study ids
  tmp <- data.frame(
    study_id = opt$study_id,
    phenotype_id = ifelse(opt$phenotype_id=="full", NA, opt$phenotype_id),
    credible_set = paste0(x %>% filter(is_cs==TRUE) %>% pull(snp), collapse=","),
    #### Nextflow working directory "work" hard coded - KEEP in mind!! ####
    path_rds = paste0(gsub("(.*)/work/.*", "\\1", getwd()), "/results/finemap/", sp_file_name, "_finemap.rds"),
    path_ind_snps = paste0(gsub("(.*)/work/.*", "\\1", getwd()), "/results/gwas_and_loci_tables/", opt$study_id, "_final_ind_snps_table.tsv")
  )
  fwrite(tmp, paste0(sp_file_name, "_coloc_info_table.tsv"),
         sep="\t", quote=F, col.names = F, na=NA)
})

cat("done!\n\n")
cat("Run-COJO rule is finished!\n")

suppressMessages(library(optparse))

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--pipeline_path", default=NULL, help="Path where Rscript lives"),
  make_option("--chr", default=NULL, help="Locus chromosome"),
  make_option("--start", default=NULL, help="Locus starting position"),
  make_option("--end", default=NULL, help="Locus ending position"),
  make_option("--phenotype_id", default=NULL, help="Trait for which the locus boundaries have been identified"),
  make_option("--dataset_gwas", default=NULL, help="GENOME-WIDE munged and aligned dataset file"),
  make_option("--mapping", default=NULL, help="Mapping file containing variants IDs matching with genotype bfile"),
  make_option("--bfile", default=NULL, help="Path and prefix name of custom LD bfiles (PLINK format .bed .bim .fam)"),
  make_option("--plink2_bin", default="/ssu/gassu/software/plink/2.00_20211217/plink2", help="Path to plink2 software"),
  make_option("--gcta_bin", default="/ssu/gassu/software/GCTA/1.94.0beta/gcta64", help="Path to GCTA software"),
  make_option("--hole", default=250000, help="Minimum pair-base distance between SNPs in different loci"),
  make_option("--cs_thresh", default=NULL, help="Percentage of credible set"),
  make_option("--outdir", default=NULL, help="Output directory"),
  make_option("--plink2_mem", default=NULL, help="Amount of RAM necessary for genotype extraction"),
  make_option("--plink2_threads", default=NULL, help="Number of threads for genotype extraction"),
  make_option("--p_label",   default=NULL, help="Label of P column"),
  make_option("--chr_label", default=NULL, help="Label of CHR column"),
  make_option("--pos_label", default=NULL, help="Label of POS column"),
  make_option("--snpid_label", default=NULL, help="Label of SNPid column"),
  make_option("--ea_label", default=NULL, help="Label of effect/minor allele column"),
  make_option("--oa_label", default=NULL, help="Label of other/non-effect allele column"),
  make_option("--eaf_label", default=NULL, help="Label of effect/minor AF column"),
  make_option("--se_label", default=NULL, help="Label of SE column"),
  make_option("--beta_label", default=NULL, help="Label of beta column"),
  make_option("--n_label", default=NULL, help="Label of sample size column"),
  make_option("--p_cojo", default=1e-04, help="P-value significant threshold for COJO (--cojo-p)"),
  make_option("--p_jumper", default=1e-04, help="P-value threshold for non-significant SNP turning significant at joint model in COJO"),
  make_option("--p_signif", default=1e-06, help="P-value significant threshold for top conditional SNP post-COJO"),
  make_option("--p_limit",  default=1e-04, help="P-value threshold for redefining locus borders post-COJO"),
  make_option("--build", default=NULL, help="Genomic build"),
  make_option("--lb_bis", default="Yes", help="Redefine locus borders post-COJO")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Source function R functions
source(paste0(opt$pipeline_path, "funs_locus_breaker_cojo_finemap_all_at_once.R"))

# parameters from config file
build <- as.character(opt$build)
redefine_region <- opt$lb_bis %in% c(TRUE, "yes", "true", "TRUE", "Yes", "1")
#Note: adding type="character", metavar="character" does not avoid changing case-sensitive Yes/No to True/False

# return input value as a string
chr.label <- sym(opt$chr_label)
pos.label <- sym(opt$pos_label)
snpid.label <- sym(opt$snpid_label)
ea.label  <- sym(opt$ea_label)
oa.label  <- sym(opt$oa_label)
eaf.label <- sym(opt$eaf_label)
se.label  <- sym(opt$se_label)
n.label   <- sym(opt$n_label)
p.label   <- sym(opt$p_label)

# converts inputs to numeric
opt$chr    <- as.numeric(opt$chr)
opt$start  <- as.numeric(opt$start)
opt$end    <- as.numeric(opt$end)

# locus name
locus_name <- paste0(opt$chr, "_", opt$start, "_", opt$end)
cat(paste("\nlocus is:", locus_name))


# condition not to run COJO on loci in HLA region
if(build == "37" & opt$chr == 6 & !(opt$end < 28477797 | opt$start > 33448354)) {
  # Create a flag file indicating an HLA signal
  flag_file <- paste0(opt$phenotype_id, "_hla_signal.txt")
  cat("HLA signal detected for", basename(opt$phenotype_id), "at this locus:", locus_name, "\n", file = flag_file)
  # stop here and move to the next locus
  #stop("The provided locus overlaps the HLA region! Stop the COJO rule here.")
  quit(status = 0)  # Exit the script silently
}


# Slightly enlarge locus by 200kb!
opt$start  <- opt$start - 100000
opt$end    <- opt$end + 100000

# to avoid killing plink job, reduce resources
opt$plink2_mem <- as.numeric(opt$plink2_mem) - 512


# reading GWAS and mapping files
dataset_gwas <- fread(opt$dataset_gwas, data.table=F)

cat(paste0("\nAdding original alleles from mapping to GWAS summary..."))

# merge map file with GWAS results
dataset_gwas <- dataset_gwas %>%
  dplyr::mutate(
    snp_map = !!snpid.label, # to report cojo results
    sdY = coloc:::sdY.est(!!se.label, !!eaf.label, !!n.label),
    type = paste0('quant') # necessary column for fine-mapping
  ) %>%
  rename(SNP = !!snpid.label)

cat(paste0("done."))

###############
# Perform cojo
###############

cat(paste0("\nRun COJO...\n\n"))
# Break dataframe in list of single rows
conditional.dataset <- cojo.ht(
  D = dataset_gwas,
  chr.label = opt$chr_label,
  pos.label = opt$pos_label,
  ea.label = opt$ea_label,
  oa.label = opt$oa_label,
  eaf.label = opt$eaf_label,
  se.label = opt$se_label,
  beta.label = opt$beta_label,
  n.label = opt$n_label,
  p.label = opt$p_label,
  locus_chr = opt$chr,
  locus_start = opt$start,
  locus_end = opt$end,
  p.cojo = as.numeric(opt$p_cojo),
  p.jumper = as.numeric(opt$p_jumper),
  bfile = opt$bfile,
  gcta.bin = opt$gcta_bin,
  plink.bin = opt$plink2_bin,
  plink.mem = opt$plink2_mem,
  plink.threads = opt$plink2_threads
)

# create folder to save outputs for each seqid separately
dir.create(paste0(opt$outdir), recursive = TRUE)
saveRDS(conditional.dataset, file=paste0(opt$outdir, "/conditional_data_", locus_name, ".rds"))


####################
# Locus breaker BIS
###################


### Repeat only on dataset that have been conditioned!!
if (redefine_region) {
  conditional.dataset$results <- lapply(conditional.dataset$results, function(x){

    ### Check if there's any SNP at p-value lower than the set threshold. Otherwise stop here
    if(isTRUE(any(x %>% pull(mlog10pC) > -log10(opt$p_signif)))){
      new_bounds <- locus.breaker(
        x,
        p.sig   = as.numeric(-log10(opt$p_signif)),
        p.limit = as.numeric(-log10(opt$p_limit)),
        hole.size = opt$hole,
        p.label   = "mlog10pC",
        chr.label = "Chr",
        pos.label = "bp")
  
      # Slightly enlarge locus by 200kb!
      new_bounds <- new_bounds %>% dplyr::mutate(start=as.numeric(start)-100000, end=as.numeric(end)+100000)

      # Remove SNPs not included in loci boundaries
      x %>% filter(bp >= new_bounds$start & bp <= new_bounds$end)
    }
  })

  ## Remove eventually empty dataframes (caused by p_thresh4 filter)
  conditional.dataset$results <- conditional.dataset$results %>% discard(is.null)
  
  ## Remove independent SNPs whose conditional data not falling in loci boundaries
  conditional.dataset$ind.snps <- conditional.dataset$ind.snps %>% filter(SNP %in% names(conditional.dataset$results))
} else {
  cat("Skipping region re-definition using locus breaker BIS function.")
}


# Quit if no signal remain after removing SNPs not included in loci boundaries
if (nrow(conditional.dataset$ind.snps) == 0) {
  flag_file <- paste0(opt$phenotype_id, "_signal_discarded.txt")
  cat("No signal remain after filter conditional datasets of", basename(opt$phenotype_id), "at locus", locus_name, "for conditional P-value <", opt$p_thresh4, "\n", file = flag_file)
  quit(status = 0)  # Exit the script silently
}

# Save updated conditional data
saveRDS(conditional.dataset, file=paste0(opt$outdir, "/conditional_data_", locus_name, "_up.rds"))

################
# Regional Plots
################

cat(paste0("done.\nTime to draw regional association plot..."))

# Plot conditioned GWAS sum stats
### have the original loci boundaries in the name, or the slightly enlarged ones?
pdf(paste0(opt$outdir, "/locus_chr", locus_name, "_conditioned_loci.pdf"), height=5*nrow(conditional.dataset$ind.snps), width=10) 
plot.cojo.ht(conditional.dataset) + patchwork::plot_annotation(paste("Locus chr", locus_name))
dev.off()


plt_loci <- plot.cojo.ht(conditional.dataset) + patchwork::plot_annotation(paste0("Locus chr", locus_name))
# Setting try() to avoid error: 
# One or both dimensions exceed the maximum (50000px). Use `options(ragg.max_dim = ...)` to change the max; Warning: May cause the R session to crash
try(suppressMessages(suppressWarnings(
  ggsave(plt_loci, filename = paste0(opt$outdir, "/locus_chr", locus_name, "_conditioned_loci.png"), height=4+2*nrow(conditional.dataset$ind.snps), width=12, dpi = 300, units = "in", limitsize = FALSE)
)))

cat("created!\n")


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
# Organize list of what needs to be saved
#########################################

cat("\nSaving independent signals...")
## Save independent association signals
fwrite(conditional.dataset$ind.snps, paste0(opt$phenotype_id, "_locus_chr", locus_name,"_ind_snps.tsv"), sep="\t", quote=F, na=NA)

cat("done.\nSave other lABF results...")

# create folder to save outputs for each seqid separately
dir.create(paste0(opt$outdir, "/finemaping/"), recursive = TRUE)

## Save lABF of each conditional dataset
lapply(finemap.res, function(x){
  sp_file_name <- paste0(opt$phenotype_id, "/finemaping/", unique(x$cojo_snp), "_locus_chr", locus_name)
  # .rds object collecting 1) lABF, 2) beta, 3) pos for all SNPs, 3) list of SNPs in the credible set
  saveRDS(x, file=paste0(sp_file_name, "_finemap.rds")) ### cojo_snp reported in the file name   #x %>% select(-cojo_snp)
  # .tsv with 1) study id and trait (if molQTL) locus info, 2) list of SNPs in the 99% credible set, 3) path and name of correspondent .rds file and 4) path and name of correspondent ind_snps.tsv table
  #  --> append each row to a master table collecting all info from processed sum stats
  ### Idea: create guidelines for generating study ids
  tmp <- data.frame(
    phenotype_id = opt$phenotype_id,
    credible_set = paste0(x %>% filter(is_cs==TRUE) %>% pull(snp), collapse=","),
    #### Nextflow working directory "work" hard coded - KEEP in mind!! ####
    path_rds      = paste0(sp_file_name, "_finemap.rds"),
    path_ind_snps = paste0(sp_file_name, "_final_ind_snps_table.tsv")
  )
  fwrite(tmp, paste0(sp_file_name, "_coloc_info_table.tsv"), sep="\t", quote=F, col.names = F, na=NA)
})

cat("done!\n\n")
cat("Run-COJO rule is finished!\n")

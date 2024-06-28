suppressMessages(library(optparse))

option_list <- list(
  make_option("--pipeline_path", default=NULL, help="Path where Rscript lives"),
  make_option("--input", default=NULL, help="Path and file name of GWAS summary statistics"),
  make_option("--p_thresh1", default=5e-08, help="Significant p-value threshold for top hits"),
  make_option("--p_thresh2", default=1e-05, help="P-value threshold for loci borders"),
  make_option("--hole", default=250000, help="Minimum pair-base distance between SNPs in different loci"),
  make_option("--study_id", default=NULL, help="Id of the study"),
  make_option("--outdir", default=NULL, help="Output directory"),
  make_option("--p_label", default=NULL, help="Label of P column"),
  make_option("--chr_label", default=NULL, help="Label of CHR column"),
  make_option("--pos_label", default=NULL, help="Label of POS column")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Source function R functions
source(paste0(opt$pipeline_path, "funs_locus_breaker_cojo_finemap_all_at_once.R"))
cat("\n\nChecking input file...")

## Throw error message - GWAS summary statistics file MUST be provided!
if(is.null(opt$input)){
  print_help(opt_parser)
  stop("Please specify the path and file name of your GWAS summary statistics in --path option", call.=FALSE)
}

cat("done!\n")

##################################
# Load in and munge GWAS sum stats
##################################

cat("\n\nReading GWAS file...")
gwas <- data.table::fread(opt$input, data.table = F)
cat("done!\n")

################
# Locus breaker
################

cat(paste0("\n\nDefining locus..."))

# function to find the index variants at each locus
check_signif <- function(x){
  ### Check if there's any SNP at p-value lower than the set threshold. Otherwise stop here
  if(any(x %>% pull(opt$p_label) > -log10(opt$p_thresh1))){
  ### Loci identification
  locus.breaker(
    x,
    p.sig     = -log10(opt$p_thresh1),
    p.limit   = -log10(opt$p_thresh2),
    hole.size = opt$hole,
    p.label   = opt$p_label,
    chr.label = opt$chr_label,
    pos.label = opt$pos_label
  )
  }else{
    stop("There is no significant signal in the current GWAS file.")
  }
} %>% discard(is.null)

cat(paste0("apply the function..."))

#list of index variants
loci_list <- check_signif(gwas)

cat(paste0("done!\n\nSaving index variants..."))

### Add study ID to the loci table. Save
loci_list$study_id <- opt$study_id

fwrite(loci_list, paste0(opt$outdir, "_loci.csv"), sep=",", quote=F, na=NA)
cat(paste0("done!\n"))

cat(paste0("\n", nrow(loci_list), " significant loci identified for ", opt$study_id, "\n"))
cat(paste0("\n\nLocus breaker is finished!\n\n"))

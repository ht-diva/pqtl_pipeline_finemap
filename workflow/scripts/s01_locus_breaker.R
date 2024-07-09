suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))


option_list <- list(
  make_option("--pipeline_path", default=NULL, help="Path where Rscript lives"),
  make_option("--input", default=NULL, help="Path and file name of GWAS summary statistics"),
  make_option("--p_thresh1", default=5e-08, help="Significant p-value threshold for top hits"),
  make_option("--p_thresh2", default=1e-05, help="P-value threshold for loci borders"),
  make_option("--hole", default=250000, help="Minimum pair-base distance between SNPs in different loci"),
  make_option("--phenotype_id", default=NULL, help="Trait for which the locus boundaries have been identified"),
  make_option("--outdir", default=NULL, help="Output directory"),
  make_option("--p_label", default=NULL, help="Label of P column"),
  make_option("--chr_label", default=NULL, help="Label of CHR column"),
  make_option("--pos_label", default=NULL, help="Label of POS column")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Source function for locus.breaker: works with -log10p
locus.breaker <- function(
    res,
    p.sig     = -log10(5e-08),
    p.limit   = -log10(1e-06),
    hole.size = 250000,
    p.label   = "LOG10P",
    chr.label = "CHROM",
    pos.label = "GENPOS"){

  res <- as.data.frame(res)
  res = res[order(as.numeric(res[, chr.label]), as.numeric(res[,pos.label])), ]
  res = res[which(res[, p.label] > p.limit), ]
  trait.res = c()

  for(j in unique(res[,chr.label])) {
    res.chr = res[which(res[, chr.label] == j), ]
    if (nrow(res.chr) > 1) {
      holes = res.chr[, pos.label][-1] - res.chr[, pos.label][-length(res.chr[,pos.label])]
      gaps = which(holes > hole.size)
      if (length(gaps) > 0) {
        for (k in 1:(length(gaps) + 1)) {
          if (k == 1) {
            res.loc = res.chr[1:(gaps[k]), ]
          }
          else if (k == (length(gaps) + 1)) {
            res.loc = res.chr[(gaps[k - 1] + 1):nrow(res.chr),
            ]
          } else {
            res.loc = res.chr[(gaps[k - 1] + 1):(gaps[k]),
            ]
          }
          if (max(res.loc[, p.label]) > p.sig) {
            start.pos = min(res.loc[, pos.label], na.rm = T)
            end.pos = max(res.loc[, pos.label], na.rm = T)
            chr = j
            best.snp = res.loc[which.max(res.loc[, p.label]),
            ]
            line.res = c(chr, start.pos, end.pos, unlist(best.snp))
            trait.res = rbind(trait.res, line.res)
          }
        }
      } else {
        res.loc = res.chr
        if (max(res.loc[, p.label]) > p.sig) {
          start.pos = min(res.loc[, pos.label], na.rm = T)
          end.pos = max(res.loc[, pos.label], na.rm = T)
          chr = j
          best.snp = res.loc[which.max(res.loc[, p.label]),
          ]
          line.res = c(chr, start.pos, end.pos, unlist(best.snp))
          trait.res = rbind(trait.res, line.res)
        }
      }
    }
    else if (nrow(res.chr) == 1) {
      res.loc = res.chr
      if (max(res.loc[, p.label]) > p.sig) {
        start.pos = min(res.loc[, pos.label], na.rm = T)
        end.pos = max(res.loc[, pos.label], na.rm = T)
        chr = j
        best.snp = res.loc[which.max(res.loc[, p.label]),
        ]
        line.res = c(chr, start.pos, end.pos, unlist(best.snp))
        trait.res = rbind(trait.res, line.res)
      }
    }
  }
  if(!is.null(trait.res)){
    trait.res = as.data.frame(trait.res, stringsAsFactors = FALSE)
    trait.res = trait.res[, -(which(names(trait.res) == chr.label))]
    names(trait.res)[1:3] = c("chr", "start", "end")
    rownames(trait.res) <- NULL
  }
  return(trait.res)
}


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

# Define a function to create a flag file indicating completion
create_flag_file <- function(outdir) {
  flag_file <- paste0(outdir, "_completed.flag")
  file.create(flag_file)
}

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
    # Create a flag file waving that there is no significant signal in the current GWAS file.
    flag_file <- paste0(opt$phenotype_id, "_no_signal.txt")
    cat("No signal detected in the GWAS of", basename(opt$phenotype_id), "at the significance level of:", opt$p_thresh1, "\n", file = flag_file)
    return(NA)  # continue to run the rest
  }
} %>% discard(is.null)

cat(paste0("apply the function..."))

#list of index variants
loci_list <- check_signif(gwas)

cat(paste0("done!\n\nSaving index variants..."))

### Add study ID to the loci table. Save
loci_list$phenotype_id <- opt$phenotype_id
fwrite(loci_list, paste0(opt$outdir, "_loci.csv"), sep=",", quote=F, na=NA)

cat(paste0("done!\n"))

cat(paste0("\n", nrow(loci_list), " significant loci identified for ", opt$phenotype_id, "\n"))
cat(paste0("\n\nLocus breaker is finished!\n\n"))


suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))

# get arguments from sbatch
option_list <- list(
  make_option("--coloc", default=NULL, help="Coloc results combined by chromosome"),
  make_option("--ld",    default=NULL, help="Pairs of coloc variants with LD r2"),
  make_option("--ofile", default=NULL, help="Output filename")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



# read inputs
ld <- read.delim(opt$ld)
res_coloc <- fread(opt$coloc)

# filter for significant coloc pairs
coloc <- res_coloc[res_coloc$PP.H4.abf > 0.8, ]
rows_with_semicolon <- coloc[grepl(";", top_cond_a)]

# Ensure data.tables
setDT(coloc)
setDT(ld)

# ID for each row
coloc[, row_id := .I]

# Funzione per espandere una riga
expand_row <- function(row) {
  
  snps_a <- unlist(strsplit(row$top_cond_a, ";\\s*"))
  snps_b <- unlist(strsplit(row$top_cond_b, ";\\s*"))
  
  # Expand all possible combinations of SNP_A and SNP_B
  expanded <- CJ(SNP_A = snps_a, SNP_B = snps_b, unique = TRUE)
  expanded[, row_id := row$row_id]
  
  return(expanded)
}

# Applichiamo su tutte le righe
coloc_expanded <- coloc[, expand_row(.SD), by = row_id]

# Mergiamo con LD
merged <- merge(
  coloc_expanded,
  ld,
  by = c("SNP_A", "SNP_B"),
  all.x = TRUE
  )

# compute r2 for each original row
r2_stats <- merged[, .(
  min_r2 = min(r2, na.rm = TRUE),
  max_r2 = max(r2, na.rm = TRUE),
  mean_r2 = mean(r2, na.rm = TRUE)
), by = row_id]

# recombine with the original coloc
coloc_ld <- merge(coloc, r2_stats, by = "row_id")

# store output
fwrite(coloc_ld, file = opt$ofile, sep = "\t", quote=F, row.names = F)

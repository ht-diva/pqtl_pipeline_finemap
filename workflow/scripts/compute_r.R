cat("\nLoading packages...\n")
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

# taking variants file as input
args <- commandArgs(trailingOnly = TRUE)

path_a <- args[1]
path_b <- args[2]
ofile  <- args[3]


path_a <- as.character(path_a)
path_b <- as.character(path_b)
ofile  <- as.character(ofile)

cat("\nLoading GWAS...\n")

pwas_a <- data.table::fread(path_a) %>% dplyr::select(SNPID, BETA)
pwas_b <- data.table::fread(path_b) %>% dplyr::select(SNPID, BETA)

cat("\nMerging packages...\n")

pwas_ab <- inner_join(
    pwas_a,
    pwas_b,
    join_by(SNPID),
    suffix = c("_a", "_b")
)

# save seqid and number of variants
seqid <- path_a %>% str_extract("seq.\\d+.\\d+")
n_a <- nrow(pwas_a)
n_b <- nrow(pwas_b)
n_overlap <- nrow(pwas_ab)

r <- cor(
    pwas_ab$BETA_a,
    pwas_ab$BETA_b,
    method = "pearson"
)

res <- data.frame(
    "seqid" = seqid,
    "nsnps_a" = n_a,
    "nsnps_b" = n_b,
    "nsnps" = n_overlap,
    "r_betas" = r, 
    "path_a" = path_a, 
    "path_b" = path_b
    )

print(res)

# store correlation
data.table::fwrite(res, file = ofile, quote = F, sep = "\t")

cat("\nCorrelation between BETA's stored here:", ofile)

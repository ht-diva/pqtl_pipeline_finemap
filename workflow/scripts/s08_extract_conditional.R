
# Required inputs to extract conditional data
# 1. identifier of protein sequence
# 2. identifier of target SNP

# For more information, please refer to the issue #26 of project page on github.

seq <- "seq.12603.87" #"seq.3054.3"
target <- "3:56849749:C:T" #"16:72114002:C:T" #"19:54759361:C:T"

# path to final(?) COJO analysis results
base <- "/scratch/dariush.ghasemi/projects/pqtl_pipeline_finemap/results/meta_filtered/"
lb_file <- "break/collected_loci_excluding_mhc.csv"
path_exm <- "seq.7128.9/conditional_data_10_115706745_116171638.rds"
path_mvp <- "/group/diangelantonio/users/giulia.pontali/file_to_share/unconditional_IVs_from_MVP.txt"


#----------------------------------------#
#-----       Prepare MVP file       -----
#----------------------------------------#

# read MVP coloc results
mvp <- data.table::fread(path_mvp) %>% head(100) %>%
  dplyr::mutate(
    seqid = str_remove_all(CHRPOS_ID, ".*_"),     # extract seqid
    locus = str_remove(locus_START_END_37, "chr") # remove 'chr' prefix from locus
    ) %>%
  dplyr::distinct(CHRPOS_ID, .keep_all = T) %>% # remove multiple associations for a target SNP 
  dplyr::select(seqid, locus, SNP, TISSUE, GENE_NAME , UNIPROT, PROTEIN_NAME) %>% #, PROTEIN_LONG_NAME
  dplyr::mutate(
    file = paste0(base, "cojo/", seqid, "/conditional_data_", locus, ".rds") # path to cojo-cond results
    )


# read final(?) merged locus breaker results
# lb_data <- read.csv(paste0(base, lb_file))

# extract locus through matching seqid and SNPID
# locus_data <- lb_data %>%
#   mutate(locus = paste(chr, start, end, sep = "_")) %>% 
#   dplyr::filter(phenotype_id %in% seq, SNPID %in% target) %>%
#   pull(locus)

# read cond data of matched locus
# cond_file <- list.files(
#   pattern = paste0("conditional_data_", locus_data, ".rds"), # pattern of output files
#   path = paste0(base, "cojo/", seq), # the path where the conditional model file seats
#   recursive = TRUE, # to show the files in subdirectories or subfolders
#   full.names = TRUE # to show full path
# )

#----------------------------------------#
#-----       Extract RDS file       -----
#----------------------------------------#

# function to:
   # 1. read RDS from path
   # 2. extract conditional data for given target SNP
   # 3. append MVP-desired columns to conditional data
grep_conditional <- function(cond_file, target){
  
  # read conditional data of corresponding locus
  cond_data <- readRDS(cond_file)
  
  # Access conditional data by target SNP identifier
  cond_data_target <- cond_data$results[[target]]
  
  mvpt <- mvp %>% dplyr::filter(SNP == target, file == cond_file)
  
  # remove unwanted columns by MVP
  cond_data_target %>% 
    #as_tibble() %>%
    dplyr::select(
      SNP, bp, refA, freq, b, se, n # freq_geno, -> to be checked with Claudia which to take
      ) %>%
  dplyr::mutate(
    seqid = mvpt$seqid,
    locus = mvpt$locus,
    TISSUE = mvpt$TISSUE,
    GENE_NAME = mvpt$GENE_NAME,
    UNIPROT = mvpt$UNIPROT,
    PROTEIN_NAME = mvpt$PROTEIN_NAME
    )
}


map2(mvp$file, mvp$SNP, grep_conditional)

# save conditional data file
write.csv(
  cond_data_target, 
  paste("conditional_data", seq, locus_data, str_replace_all(target,":","_"), ".csv", sep = "_"),
  row.names = F, quote = F
  )








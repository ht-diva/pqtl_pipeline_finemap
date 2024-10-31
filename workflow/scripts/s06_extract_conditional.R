
# required inputs to extract conditional data
# 1. identifier of protein sequence
# 2. identifier of target SNP

seq <- "seq.12603.87" #"seq.3054.3"
target <- "3:56849749:C:T" #"16:72114002:C:T" #"19:54759361:C:T"

# path to final(?) COJO analysis results
base <- "/scratch/dariush.ghasemi/projects/pqtl_pipeline_finemap/results/meta_filtered/"
lb_file <- "break/collected_loci_excluding_mhc.csv"

# read final(?) merged locus breaker results
lb_data <- read.csv(paste0(base, lb_file))

# extract locus through matching seqid and SNPID
locus_data <- lb_data %>%
  mutate(locus = paste(chr, start, end, sep = "_")) %>% 
  dplyr::filter(phenotype_id %in% seq, SNPID %in% target) %>%
  pull(locus)

# read cond data of matched locus
cond_file <- list.files(
  pattern = paste0("conditional_data_", locus_data, ".rds"), # pattern of output files
  path = paste0(base, "cojo/", seq), # the path where the conditional model file seats
  recursive = TRUE, # to show the files in subdirectories or subfolders
  full.names = TRUE # to show full path
)

# read conditional data of corresponding locus
cond_data <- readRDS(cond_file)

# Access conditional data by target SNP identifier
cond_data_target <- cond_data$results[[target]]

# save conditional data file
write.csv(
  cond_data_target, 
  paste("conditional_data", seq, locus_data, str_replace_all(target,":","_"), ".csv", sep = "_"),
  row.names = F, quote = F
  )


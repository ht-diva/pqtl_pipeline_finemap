suppressMessages(library(optparse))

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--pipeline_path", default=NULL, help="Path where Rscript lives"),
  make_option("--input", default=NULL, help="Path and filename of master coloc table produced by individual traits preprocessing")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Source function R functions
source(paste0(opt$pipeline_path, "funs_locus_breaker_cojo_finemap_all_at_once.R"))

# List all cs variants

###### TO REMOVE
tb <- rbind(
  fread("results/final/ATAC_chr22_coloc_info_master_table.tsv", data.table=F),
  fread("results/final/Cardioembolic_Stroke_Eur_Mishra_2022_Nature_coloc_info_master_table.tsv", data.table=F))
########
tb <- fread(opt$input, data.table=F)

tb <- tb %>% mutate(cs_name=paste0("cs", seq(1,nrow(tb))))
cs_list <- lapply(tb$credible_set, function(x){ unlist(str_split(x, ",")) })
names(cs_list) <- tb$cs_name

# Get unique elements present in all text vectors
all_elements <- unique(unlist(cs_list))

# Create a sparse matrix with zeros
matrix_table_sparse <- Matrix(
  0,
  nrow = length(cs_list),
  ncol = length(all_elements),
  sparse = TRUE,
  dimnames = list(names(cs_list), all_elements)
)

# Fill in the matrix with 1 where the text vector contains the element
for (i in seq_along(cs_list)) {
  matrix_table_sparse[i, cs_list[[i]]] <- 1
}


# Multiply the transpose of the matrix with itself
matrix_product <- matrix_table_sparse %*% t(matrix_table_sparse)

# Set the diagonal to zero (as each vector will have 1s on diagonal)
diag(matrix_product) <- 0

shared_elements <- as.data.frame(which(matrix_product > 0, arr.ind = TRUE))

# For each unique SNP, list of cs having it
shared_elements <- as.data.frame(shared_elements) %>%
  mutate(row=rownames(shared_elements)) %>%
  rowwise() %>%
  mutate(col=paste0(c("cs", col), collapse="")) %>%
  mutate(check=paste0(sort(c(row, col)), collapse="")) %>%
  distinct(check, .keep_all = T) %>%
  select(-check)

### Credible sets not overlapping with anything -- TO KEEP FOR LATER
# cs_alone <- unique(c(setdiff(names(cs_list), shared_elements$row), setdiff(names(cs_list), shared_elements$col)))

# Retrieve trait info
coloc_combo <- shared_elements %>%
  rename(t1=row, t2=col) %>%
  left_join(tb %>% select(-credible_set), by=c("t1"="cs_name")) %>%
  dplyr::rename(t1_study_id=study_id, t1_phenotype_id=phenotype_id, t1_path=path) %>%
#  mutate(
#    t1_locus=gsub(".*/.*_locus_(chr\\d+_\\d+_\\d+).*", "\\1", t1_path),
#    t1_top_snp=gsub(".*/.*_top_snp_(.*)_finemap.rds", "\\1", t1_path)
#  ) %>%
  left_join(tb %>% select(-credible_set), by=c("t2"="cs_name")) %>%
  dplyr::rename(t2_study_id=study_id, t2_phenotype_id=phenotype_id, t2_path=path) #%>%
#  mutate(
#    t2_locus=gsub(".*/.*_locus_(chr\\d+_\\d+_\\d+).*", "\\1", t2_path),
#    t2_top_snp=gsub(".*/.*_top_snp_(.*)_finemap.rds", "\\1", t2_path)
#  )


# Remove pair testing different conditional dataset for the same trait (study_id + phenotype_id)
coloc_combo <- coloc_combo %>%
  filter(t1_study_id != t2_study_id | t1_study_id==t2_study_id & t1_phenotype_id != t2_phenotype_id) %>%
  select(-t1, -t2)

fwrite(coloc_combo, "coloc_pairwise_guide_table.tsv", quote=F, na=NA, sep="\t")

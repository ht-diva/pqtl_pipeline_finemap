suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(Matrix))

# Get arguments specified in the sbatch
option_list <- list(
  #make_option("--pipeline_path", default=NULL, help="Path where Rscript lives"),
  make_option("--coloc_info_table", default=NULL, help="Path and filename of master coloc table produced by individual traits pre-processing"),
#  make_option("--coloc_id", default=NULL, help="Id code to univocally identify the colocalisation analysis"),
  make_option("--chr_cs", default=NULL, help="Chromosomes to split the analysis in"),
#  make_option("--chunk_size", default=NULL, help="Row number of each chunk the output file will be broken into")
  make_option("--ofile", default=NULL, help="Output filename")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# List all cs variants
tb <- data.table::fread(opt$coloc_info_table, header = T) %>% dplyr::filter(chr == opt$chr_cs)

nrows_tb <- nrow(tb)

# Filter by chr
tb <- tb  %>% dplyr::mutate(cs_name=paste0("cs", seq(1,nrows_tb)))

cs_list <- strsplit(tb$credible_set, ",", fixed=TRUE) # If TRUE match split exactly, otherwise use regular expressions - hoping is more memory efficient! https://stackoverflow.com/questions/55919893/more-memory-efficient-way-than-strsplit-to-split-a-string-into-two-in-r
names(cs_list) <- tb$cs_name

# Get unique elements present in all text vectors
all_elements <- unique(unlist(cs_list))

# Check that there's at list a case of overlapping! Otherwise no point in performing all this - FOR THE MOMENT. Still useful to add it to the massive matrix? But will decide later

if(length(unlist(cs_list)) > length(all_elements)){
  
  element_indices <- setNames(seq_along(all_elements), all_elements)
  
  # Create a sparse matrix with zeros
  matrix_sparse <- Matrix::Matrix(
    0,
    nrow = length(cs_list),
    ncol = length(all_elements),
    sparse = TRUE,
    dimnames = list(names(cs_list), all_elements)
  )
  
  
  # Create a matrix of indices for each row
  row_indices <- lapply(cs_list, function(vec) {
    element_indices[vec]
  })
  
  
  matrix_sparse[cbind(rep(seq_along(cs_list), lengths(row_indices)), unlist(row_indices))] <- 1
  
  # Fill in the matrix with 1 where the text vector contains the element
  
  # Multiply the transpose of the matrix with itself
  start <- Sys.time()
  matrix_product <- matrix_sparse %*% Matrix::t(matrix_sparse)

  
  # Set the diagonal to zero (as each vector will have 1s on diagonal)
  diag(matrix_product) <- 0
  
  shared_elements <- Matrix::which(matrix_product > 0, arr.ind = TRUE)
  shared_elements_df <- as.data.frame(shared_elements)
  
  # Map row and column indices to credible set names
  shared_elements_df$row <- rownames(shared_elements)
  shared_elements_df$col <- colnames(matrix_product)[shared_elements_df$col]
  
  # Keep unique pairs
  shared_elements_unique <- unique(Matrix::t(apply(shared_elements_df, 1, sort)))
  colnames(shared_elements_unique) <- c("t1", "t2")
  
  # Retrieve trait info
  coloc_combo <- merge(
    shared_elements_unique,
    tb %>% dplyr::select(-credible_set),
    by.x = "t1",
    by.y = "cs_name"
  )
  
  coloc_combo <- merge(
    coloc_combo,
    tb %>% dplyr::select(-credible_set),
    by.x = "t2",
    by.y = "cs_name",
    suffixes = c("_t1", "_t2")
  )
  
  # Remove pair testing different conditional dataset for the same trait (study_id + phenotype_id)
  coloc_combo <- coloc_combo %>%
    dplyr::filter(phenotype_id_t1 != phenotype_id_t2) %>%
    dplyr::select(-t1, -t2, -chr_t1, -chr_t2) %>%
    dplyr::rename(
      t1_phenotype_id = phenotype_id_t1,
      t2_phenotype_id = phenotype_id_t2,
      t1_path_rds = path_rds_t1,
      t2_path_rds = path_rds_t2,
      t1_path_ind_snps = path_ind_snps_t1,
      t2_path_ind_snps = path_ind_snps_t2
      )
  
  end <- Sys.time()
  cat("\nAnalysis time: ",end-start)
  
  cat("\nsaving results to: ", paste0(opt$ofile))
  
  # Save the table as input of next rule
  data.table::fwrite(coloc_combo, file = paste0(opt$ofile), quote=F, na=NA, sep="\t")
}

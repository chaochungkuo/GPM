
"Module for converting Seurat objects stored as 'RDS' objects to Scanpy-readable h5ad. The created file has same name and saved in the same directory as the input file."

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")


    # Define the packages to check and install if necessary
packages <- c("GenomeInfoDb", "GenomeInfoDbData")

# Identify packages that are not already installed
packages_to_install <- packages[!(packages %in% installed.packages()[, "Package"])]

# Install the missing packages without prompting any questions and suppress messages and warnings
if (length(packages_to_install) > 0) {
    suppressMessages(suppressWarnings(BiocManager::install(packages_to_install, ask = FALSE)))
} else {
    message("All specified packages are already installed.")
}


suppressMessages(library(SingleCellExperiment, quietly = T, verbose = F, warn.conflicts = F))
suppressMessages(library(zellkonverter, quietly = T, verbose = F, warn.conflicts = F))
suppressMessages(library(optparse, quietly = T, verbose = F, warn.conflicts = F))


convert_h5ad_to_rds <- function(h5ad_file) {
  sce_file <- sub("\\.h5ad$", ".rds", h5ad_file)
  cat(sprintf("Converting %s to %s...\n", h5ad_file, sce_file))
  
  sce <- readH5AD(h5ad_file)
  saveRDS(sce, sce_file, compress = FALSE)
  
  cat(sprintf("Conversion complete: %s\n", sce_file))
}

convert_rds_to_h5ad <- function(rds_file) {
  h5ad_file <- sub("\\.rds$", ".h5ad", rds_file)
  cat(sprintf("Converting %s to %s...\n", rds_file, h5ad_file))
  
  sce <- readRDS(rds_file)
  writeH5AD(sce, h5ad_file, overwrite = TRUE)
  
  cat(sprintf("Conversion complete: %s\n", h5ad_file))
}

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input file (.h5ad or .rds)", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file must be specified.", call. = FALSE)
}

input_file <- opt$input

if (grepl("\\.h5ad$", input_file)) {
  convert_h5ad_to_rds(input_file)
} else if (grepl("\\.rds$", input_file)) {
  convert_rds_to_h5ad(input_file)
} else {
  stop("Unsupported file type. Please provide a .h5ad or .rds file.", call. = FALSE)
}

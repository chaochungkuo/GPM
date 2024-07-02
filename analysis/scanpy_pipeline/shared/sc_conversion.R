

# Set the download URL
dn_url <- "https://bioconductor.org/packages/3.18/data/annotation/src/contrib/GenomeInfoDbData_1.2.11.tar.gz"

print(dn_url)

# Create a temporary directory (replace with your preferred method)
tmp_dir <- tempdir()  # This creates a temporary directory


# Construct the download path
dn_path <- file.path(tmp_dir, "GenomeInfoDbData_1.2.11.tar.gz")

# Download the file using download.file()
download.file(url = dn_url, destfile = dn_path, method = "wget")  # Use curl for efficiency
system(paste0("R CMD INSTALL ", dn_path))


suppressMessages(library(SingleCellExperiment, quietly = T, verbose = F, warn.conflicts = F))
suppressMessages(library(zellkonverter, quietly = T, verbose = F, warn.conflicts = F))
suppressMessages(library(optparse, quietly = T, verbose = F, warn.conflicts = F))


convert_h5ad_to_rds <- function(h5ad_file, output = NULL) {
  if (is.null(output)) {
  sce_file <- sub("\\.h5ad$", ".rds", h5ad_file)
  } else {
    sce_file <- output
  }
  cat(sprintf("Converting %s to %s...\n", h5ad_file, sce_file))
  
  sce <- readH5AD(h5ad_file)
  saveRDS(sce, sce_file, compress = FALSE)
  
  cat(sprintf("Conversion complete: %s\n", sce_file))
}

convert_rds_to_h5ad <- function(rds_file, output = NULL) {

  if (is.null(output)) {
    h5ad_file <- sub("\\.rds$", ".h5ad", rds_file)
  } else {
    h5ad_file <- output
  }
  cat(sprintf("Converting %s to %s...\n", rds_file, h5ad_file))
  sce <- readRDS(rds_file)
  writeH5AD(sce, h5ad_file)
  cat(sprintf("Conversion complete: %s\n", h5ad_file))
}

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input file (.h5ad or .rds)", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output file (.h5ad or .rds)", metavar = "FILE"),
  make_option(c("-h", "--help"), action = "store_true", default = FALSE),
  

)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file must be specified.", call. = FALSE)
}

input_file <- opt$input
output_file <- opt$output

if (grepl("\\.h5ad$", input_file)) {
  convert_h5ad_to_rds(input_file, output_file)
} else if (grepl("\\.rds$", input_file)) {
  convert_rds_to_h5ad(input_file, output_file)
} else {
  stop("Unsupported file type. Please provide a .h5ad or .rds file.", call. = FALSE)
}

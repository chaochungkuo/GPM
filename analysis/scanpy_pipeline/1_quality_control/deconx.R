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

suppressMessages(library(optparse))
suppressMessages(library(decontX))
suppressMessages(library(zellkonverter))
suppressMessages(library(SingleCellExperiment))

# Define command-line options
option_list <- list(
  make_option(c("-s", "--sce"), type = "character", default = NULL, help = "Path to the SCE RDS file", metavar = "character"),
  make_option(c("-r", "--raw"), type = "character", default = NULL, help = "Path to the raw RDS file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "Output file path for the saved H5AD file", metavar = "character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$sce) || is.null(opt$raw) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Please provide the paths to the SCE file, raw file, and the output file.", call. = FALSE)
}

# Load the SCE and raw files
sce <- readH5AD(opt$sce)
raw <- readH5AD(opt$raw)

# Rename assays
assays(sce, withDimnames = FALSE) <- setNames(assays(sce), gsub("X", "counts", assayNames(sce)))
assays(raw, withDimnames = FALSE) <- setNames(assays(raw), gsub("X", "counts", assayNames(raw)))

# Apply decontX
sce <- decontX(x = sce, background = raw)

# Save the result as an H5AD file
writeH5AD(sce = sce, file = opt$output)

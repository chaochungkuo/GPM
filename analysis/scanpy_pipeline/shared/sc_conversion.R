a <- installed.packages()
packages <- a[, 1]

### ----------------------------------------------------------------------###
###                       Install required packages                      ###
### ----------------------------------------------------------------------###

if (!is.element("GenomeInfoDbData", packages)) {
  dn_url <- "https://bioconductor.org/packages/3.18/data/annotation/src/contrib/GenomeInfoDbData_1.2.11.tar.gz"
  tmp_dir <- tempdir() # This creates a temporary directory
  dn_path <- file.path(tmp_dir, "GenomeInfoDbData_1.2.11.tar.gz")
  download.file(url = dn_url, destfile = dn_path, method = "wget") # Use curl for efficiency
  system(paste0("R CMD INSTALL ", dn_path))
}

suppressMessages(library(SingleCellExperiment, quietly = T, verbose = F, warn.conflicts = F))
suppressMessages(library(zellkonverter, quietly = T, verbose = F, warn.conflicts = F))
suppressMessages(library(optparse, quietly = T, verbose = F, warn.conflicts = F))


### ----------------------------------------------------------------------###
###                         Define CLI options                           ###
### ----------------------------------------------------------------------###

option_list <- list(
  make_option(c("-i", "--input"),
    type = "character", default = NULL,
    help = "Input file (.h5ad or .rds)", metavar = "FILE"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = NULL,
    help = "Output file (.h5ad, .rds or .cloupe)", metavar = "FILE"
  ),
  make_option(c("-f", "--from"),
    type = "character", default = NULL,
    help = "Input file format ('seurat, SCE, or anndata')", metavar = "FORMAT"
  ),
  make_option(c("-t", "--to"),
    type = "character", default = NULL,
    help = "Output file format ('seurat, SCE, anndata, or loupe')", metavar = "FORMAT"
  )
)


### ----------------------------------------------------------------------###
###                         Parse & validate CLI arguments                ###
### ----------------------------------------------------------------------###



opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file must be specified.", call. = FALSE)
}

if (is.null(opt$from)) {
  print_help(opt_parser)
  stop("Input format must be specified.", call. = FALSE)
}

if (is.null(opt$to)) {
  print_help(opt_parser)
  stop("Output format must be specified.", call. = FALSE)
}



input_file <- opt$input
output_file <- opt$output
output_format <- opt$to
base_name <- basename(input_file)



### ----------------------------------------------------------------------###
###                         Conversion functions                         ###
### ----------------------------------------------------------------------###

conver_sce_to_seurat <- function(sce_obj) {
  suppressMessages(library(Seurat, quietly = T, verbose = F, warn.conflicts = F))
  seurat_obj <- Seurat::as.Seurat(sce_obj)
  return(seurat_obj)
}

# Works on linux as well
convert_seurat_to_loupeR <- function(seurat_obj) {
  # From: https://github.com/10xGenomics/loupeR

  if (!is.element("loupeR", packages)) {
    # install platform specific source package
    os <- sub("Darwin", "macOS", "Linux", Sys.info()["sysname"])
    url <- paste0("https://github.com/10XGenomics/loupeR/releases/latest/download/loupeR_", os, ".tar.gz")
    install.packages(url, repos = NULL, type = "source")
  }
  suppressMessages(library(loupeR, quietly = T, verbose = F, warn.conflicts = F))

  loupeR_obj <- loupeR::create_loupe_from_seurat(seurat_obj)
  return(loupeR_obj)
}

convert_h5ad_to_rds <- function(h5ad_file, output = NULL) {
  if (is.null(output)) {
    sce_file <- sub("\\.h5ad$", ".rds", h5ad_file)
  } else {
    sce_file <- output
  }
  cat(sprintf("Converting %s to %s...\n", h5ad_file, sce_file))

  sce <- readH5AD(h5ad_file)

  if (output_format == "SCE") {
    saveRDS(sce, sce_file, compress = FALSE)
  } else if (output_format == "seurat") {
    seurat_obj <- conver_sce_to_seurat(sce)
    saveRDS(seurat_obj, sce_file, compress = FALSE)
  } else if (output_format == "loupe") {
    seurat_obj <- conver_sce_to_seurat(sce)
    loupeR_obj <- convert_seurat_to_loupeR(seurat_obj)
  } else {
    stop("Unsupported output format. Please provide 'SingleCellExperiment', 'Seurat', or 'LoupeR'.", call. = FALSE)
  }


  cat(sprintf("Conversion complete: %s\n", sce_file))
}

### ----------------------------------------------------------------------###
###                       To h5ad Convertor functions                     ###
### ----------------------------------------------------------------------###

# TODO: Add support for Seurat to SCE conversion
convert_sce_to_h5ad <- function(rds_file, output = NULL) {
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

convert_seurat_to_h5ad <- function(rds_file, output = NULL) {
  if (is.null(output)) {
    h5ad_file <- sub("\\.rds$", ".h5ad", rds_file)
  } else {
    h5ad_file <- output
  }

  cat(sprintf("Converting %s to %s...\n", rds_file, h5ad_file))
  seurat_obj <- readRDS(rds_file)
  sce <- Seurat::as.SingleCellExperiment(seurat_obj)
  writeH5AD(sce, h5ad_file)
  cat(sprintf("Conversion complete: %s\n", h5ad_file))
}


### ----------------------------------------------------------------------###
###                         Convertor dispatch                           ###
### ----------------------------------------------------------------------###

if (opt$from == "SCE" && opt$to == "h5ad") {
  convert_sce_to_h5ad(input_file, output_file)
} else if (opt$from == "seurat" && opt$to == "h5ad") {
  convert_seurat_to_h5ad(input_file, output_file)
} else if (opt$from == "h5ad" && opt$to == "rds") {
  convert_h5ad_to_rds(input_file, output_file)
} else {
  stop("Unsupported conversion. Please provide 'SCE' to 'anndata', 'seurat' to 'anndata', or 'anndata' to 'seurat'.", call. = FALSE)
}

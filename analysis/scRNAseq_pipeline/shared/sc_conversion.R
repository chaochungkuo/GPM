a <- installed.packages()
packages <- a[, 1]

### ----------------------------------------------------------------------###
###                       Install required packages                      ###
### ----------------------------------------------------------------------###

if (!is.element("GenomeInfoDbData", packages)) {
  dn_url <- "https://mghp.osn.xsede.org/bir190004-bucket01/archive.bioconductor.org/packages/3.18/data/annotation/src/contrib/GenomeInfoDbData_1.2.11.tar.gz"
  tmp_dir <- tempdir() # This creates a temporary directory
  dn_path <- file.path(tmp_dir, "GenomeInfoDbData_1.2.11.tar.gz")
  download.file(url = dn_url, destfile = dn_path, method = "curl") # Use curl for efficiency
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
###                         Parse CLI arguments                           ###
### ----------------------------------------------------------------------###



opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)



### ----------------------------------------------------------------------###
###                         Define global variables                       ###
### ----------------------------------------------------------------------###

# # accounting for re-map of file paths in docker
# if (file.exists("/.dockerenv")) {
#   docker <- TRUE
# } else {
#   docker <- FALSE
# }


# if (docker == TRUE) {
#   input_file <- paste0("/app/host", opt$input)
#   output_file <- paste0("/app/host", opt$output)
# } else {
#   input_file <- opt$input
#   output_file <- opt$output
# }

input_file <- opt$input
output_file <- opt$output
input_format <- opt$from
output_format <- opt$to



### ----------------------------------------------------------------------###
###                         Validate CLI arguments                           ###
### ----------------------------------------------------------------------###


if (is.null(input_file)) {
  print_help(opt_parser)
  stop("Input file must be specified.", call. = FALSE)
}

if (is.null(output_file)) {
  print_help(opt_parser)
  stop("Input format must be specified.", call. = FALSE)
}

if (is.null(output_format)) {
  print_help(opt_parser)
  stop("Output format must be specified.", call. = FALSE)
}

if (is.null(input_format)) {
  print_help(opt_parser)
  stop("Input format must be specified.", call. = FALSE)
}


### ----------------------------------------------------------------------###
###                         Conversion functions                         ###
### ----------------------------------------------------------------------###


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

  loupeR_obj <- loupeR::create_loupe_from_seurat(seurat_obj, output_dir = dirname(output_file), output_name = basename(output_file))
  return(loupeR_obj)
}

convert_sce_to_seurat <- function(sce_obj) {
  suppressMessages(library(Seurat, quietly = T, verbose = F, warn.conflicts = F))
  seurat_obj <- Seurat::as.Seurat(sce_obj)
  return(seurat_obj)
}

convert_seurat_to_sce <- function(seurat_obj) {
  sce_obj <- Seurat::as.SingleCellExperiment(seurat_obj)
  return(sce_obj)
}



convert_h5ad_to_sce <- function(h5ad_file, output) {
  cat(sprintf("Converting %s to %s...\n", h5ad_file, output))

  sce <- readH5AD(h5ad_file)
  cat(sprintf("Conversion complete: %s\n", output))
  return(sce)
}

convert_sce_to_h5ad <- function(rds_file, output) {
  cat(sprintf("Converting %s to %s...\n", rds_file, output))
  sce <- readRDS(rds_file)
  writeH5AD(sce, h5ad_file, compression = "none")
  cat(sprintf("Conversion complete: %s\n", output))
}



### ----------------------------------------------------------------------###
###                         Convertor dispatch                            ###
### ----------------------------------------------------------------------###

if (input_format == "SCE") {
  if (output_format == "seurat") {
    sce_obj <- readRDS(input_file)
    seurat_obj <- convert_sce_to_seurat(sce_obj)
    saveRDS(seurat_obj, output_file, compress = FALSE)
  } else if (output_format == "anndata") {
    sce_obj <- readRDS(input_file)
    writeH5AD(sce_obj, output_file, compression = "none")
  } else if (output_format == "loupe") {
    sce_obj <- readRDS(input_file)
    convert_seurat_to_loupeR(convert_sce_to_seurat(sce_obj))
  } else {
    stop("Conversion from SCE to ", output_format, " is not supported.", call. = FALSE)
  }
}

if (input_format == "seurat") {
  if (output_format == "SCE") {
    seurat_obj <- readRDS(input_file)
    sce_obj <- convert_seurat_to_sce(seurat_obj)
    saveRDS(sce_obj, output_file, compress = FALSE)
  } else if (output_format == "anndata") {
    seurat_obj <- readRDS(input_file)
    writeH5AD(convert_seurat_to_sce(seurat_obj), output_file, compression = "none")
  } else if (output_format == "loupe") {
    seurat_obj <- readRDS(input_file)
    convert_seurat_to_loupeR(seurat_obj)
  } else {
    stop("Conversion from Seurat to ", output_format, " is not supported.", call. = FALSE)
  }
}


if (input_format == "anndata") {
  if (output_format == "SCE") {
    sce_obj <- convert_h5ad_to_sce(input_file, output_file)
    saveRDS(sce_obj, output_file, compress = FALSE)
  } else if (output_format == "seurat") {
    sce_obj <- convert_h5ad_to_sce(input_file, output_file)
    seurat_obj <- convert_sce_to_seurat(sce_obj)
    saveRDS(seurat_obj, output_file, compress = FALSE)
  } else if (output_format == "loupe") {
    sce_obj <- convert_h5ad_to_sce(input_file, output_file)
    convert_seurat_to_loupeR(convert_sce_to_seurat(sce_obj))
  } else {
    stop("Conversion from anndata to ", output_format, " is not supported.", call. = FALSE)
  }
}

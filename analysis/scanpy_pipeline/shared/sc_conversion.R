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

if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Output format must be specified.", call. = FALSE)
}



input_file <- opt$input
output_file <- opt$output
output_format <- opt$to
base_name <- basename(input_file)
dir_name <- dirname(input_file)



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

  loupeR_obj <- loupeR::create_loupe_from_seurat(seurat_obj)
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
  cat(sprintf("Converting %s to %s...\n", h5ad_file, sce_file))

  sce <- readH5AD(h5ad_file)
  cat(sprintf("Conversion complete: %s\n", sce_file))
  return(sce)
}

convert_sce_to_h5ad <- function(rds_file, output) {

  cat(sprintf("Converting %s to %s...\n", rds_file, h5ad_file))
  sce <- readRDS(rds_file)
  writeH5AD(sce, h5ad_file)
  cat(sprintf("Conversion complete: %s\n", h5ad_file))

}



### ----------------------------------------------------------------------###
###                         Convertor dispatch                           ###
### ----------------------------------------------------------------------###

if (opt$from == "SCE") {

  if (opt$to == "seurat") {
    sce_obj <- readRDS(input_file)
    seurat_obj <- convert_sce_to_seurat(sce_obj)
    saveRDS(seurat_obj, output_file)

  } else if (opt$to == "anndata") {
    sce_obj <- readRDS(input_file)
    writeH5AD(sce_obj, output_file)

  } else if (opt$to == "loupe") {
    sce_obj <- readRDS(input_file)
    loupe_obj <- convert_seurat_to_loupeR(convert_sce_to_seurat(sce_obj))
    saveRDS(loupe_obj, output_file)

  } else {
    stop("Conversion from SCE to ", opt$to, " is not supported.", call. = FALSE)
  }
}

if(opt$from == "seurat") {

  if (opt$to == "SCE") {
    seurat_obj <- readRDS(input_file)
    sce_obj <- convert_seurat_to_sce(seurat_obj)
    saveRDS(sce_obj, output_file)

  } else if (opt$to == "anndata") {
    seurat_obj <- readRDS(input_file)
    writeH5AD(convert_seurat_to_sce(seurat_obj), output_file)

  } else if (opt$to == "loupe") {
    seurat_obj <- readRDS(input_file)
    loupe_obj <- convert_seurat_to_loupeR(seurat_obj)
    saveRDS(loupe_obj, output_file)

  } else {
    stop("Conversion from Seurat to ", opt$to, " is not supported.", call. = FALSE)
  }
}


if(opt$from == "anndata") {

  if (opt$to == "SCE") {
    sce_obj <- convert_h5ad_to_sce(input_file, output_file)
    saveRDS(sce_obj, output_file)

  } else if (opt$to == "seurat") {
    sce_obj <- convert_h5ad_to_sce(input_file, output_file)
    seurat_obj <- convert_sce_to_seurat(sce_obj)
    saveRDS(seurat_obj, output_file)

  } else if (opt$to == "loupe") {
    sce_obj <- convert_h5ad_to_sce(input_file, output_file)
    loupe_obj <- convert_seurat_to_loupeR(convert_sce_to_seurat(sce_obj))
    saveRDS(loupe_obj, output_file)

  } else {
    stop("Conversion from anndata to ", opt$to, " is not supported.", call. = FALSE)
  }
}

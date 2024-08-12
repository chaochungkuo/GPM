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


### ----------------------------------------------------------------------###
###                       LoupeR                                          ###
### ----------------------------------------------------------------------###


eula_create <- function() {
  dir.create(eula_data_dir(), showWarnings = FALSE, recursive = TRUE)
  file.create(eula_lock_file())
}


#' Path to directory that holds EULA agreement lock file
#' @noRd
eula_data_dir <- function() {
  tools::R_user_dir("loupeR", "data")
}

#' Path to EULA agreement lock file
#' @noRd
eula_lock_file <- function() {
  file.path(eula_data_dir(), "eula_agreement")
}

if (!is.element("loupeR", packages)) {
  # install platform specific source package
  os <- sub("Darwin", "macOS", "Linux", Sys.info()["sysname"])
  url <- paste0("https://github.com/10XGenomics/loupeR/releases/latest/download/loupeR_", os, ".tar.gz")
  install.packages(url, repos = NULL, type = "source")
  eula_create()
  loupeR::setup()
}


### ----------------------------------------------------------------------###
###                       Basilisk Env                                    ###
### ----------------------------------------------------------------------###

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("basilisk.utils")


library(zellkonverter)
library(basilisk)
library(basilisk.utils)

basilisk.utils::installConda()
system("conda install -n base conda-libmamba-solver")
system("conda config --set solver libmamba")
env <- zellkonverter::zellkonverterAnnDataEnv()
setupBasiliskEnv(envpath = paste0(Sys.getenv("HOME"), "/.cache/R/basilisk/1.14.1/zellkonverter/1.12.1/", env@envname),
                 packages = env@packages,
                 channels = env@channels,
                 pip = env@pip, 
                 paths = env@paths
                 )

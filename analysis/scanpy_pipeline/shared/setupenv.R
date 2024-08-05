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

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("basilisk.utils")


library(zellkonverter)
library(basilisk)
library(basilisk.utils)

basilisk.utils::installConda()
env = zellkonverter::zellkonverterAnnDataEnv()
setupBasiliskEnv(envpath=paste0(Sys.getenv("HOME"),"/.cache/R/basilisk/1.14.1/zellkonverter/1.12.1/", env@envname), packages= env@packages, channels = env@channels, pip=env@pip, paths = env@paths)
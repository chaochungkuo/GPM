package_list = installed.packages()

# Function to process each package
process_package <- function(package, package_file, url) {
    # Check if package is already installed
    if (package %in% rownames(package_list)) {
        cat(package, "is already installed\n")
        return(TRUE)
    }
  if (!file.exists(package_file)) {
    cat(package_file, "not found. Downloading...\n")
    tryCatch({
      download.file(url, package_file, mode = "wb")
    }, error = function(e) {
      cat("Failed to download", package_file, "\n")
      return(FALSE)
    })
  }
  
  cat("Installing", package_file, "...\n")
  install.packages(package_file, repos = NULL, type = "source")
  
  if (!file.exists(package_file)) {
    cat("Removing downloaded", package_file, "\n")
    file.remove(package_file)
  }
  
  return(TRUE)
}

packages <- c(
  "GenomeInfoDbData",
  "GO.db",
  "HDO.db",
  "org.Hs.eg.db",
  "org.Mm.eg.db"
)

# List of packages and their URLs
packages_files <- c(
  "GenomeInfoDbData_1.2.10.tar.gz",
  "GO.db_3.17.0.tar.gz",
  "HDO.db_0.99.1.tar.gz",
  "org.Hs.eg.db_3.17.0.tar.gz",
  "org.Mm.eg.db_3.17.0.tar.gz"
)

# Read URLs from download.txt
urls <- c("https://bioconductor.org/packages/3.17/data/annotation/src/contrib/GenomeInfoDbData_1.2.10.tar.gz",
        "https://bioconductor.org/packages/3.17/data/annotation/src/contrib/GO.db_3.17.0.tar.gz",
        "https://bioconductor.org/packages/3.17/data/annotation/src/contrib/HDO.db_0.99.1.tar.gz",
        "https://bioconductor.org/packages/3.17/data/annotation/src/contrib/org.Mm.eg.db_3.17.0.tar.gz",
        "https://bioconductor.org/packages/3.17/data/annotation/src/contrib/org.Hs.eg.db_3.17.0.tar.gz"
        )

# Process each package
for (i in seq_along(packages)) {
  process_package(packages[i], packages_files[i],  urls[i])
}

if ("netdiffuseR" %in% rownames(package_list)) {
    cat("netdiffuseR", "is already installed\n")
}else {
   install.packages("netdiffuseR", repos='http://cran.us.r-project.org')

}

if ("CrossTalkeR" %in% rownames(package_list)) {
    cat("CrossTalkeR", "is already installed\n")
}else {
system("wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=115M510CCnjzwXVGd3_B3RdcntvP0UPyw' -O CrossTalkeR.tar.gz")
system("R CMD INSTALL CrossTalkeR.tar.xz")
file.remove("CrossTalkeR.tar.gz")
}
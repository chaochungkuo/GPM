

# Function to process each package
process_package <- function(package, url) {
  if (!file.exists(package)) {
    cat(package, "not found. Downloading...\n")
    tryCatch({
      download.file(url, package, mode = "wb")
    }, error = function(e) {
      cat("Failed to download", package, "\n")
      return(FALSE)
    })
  }
  
  cat("Installing", package, "...\n")
  install.packages(package, repos = NULL, type = "source")
  
  if (!file.exists(package)) {
    cat("Removing downloaded", package, "\n")
    file.remove(package)
  }
  
  return(TRUE)
}

# List of packages and their URLs
packages <- c(
  "GenomeInfoDbData_1.2.10.tar.gz",
  "GO.db_3.17.0.tar.gz",
  "HDO.db_0.99.1.tar.gz",
  "org.Hs.eg.db_3.17.0.tar.gz",
  "org.Mm.eg.db_3.17.0.tar.gz"
  
)

# Read URLs from download.txt
urls <- readLines("download.txt")

if (length(packages) != length(urls)) {
  stop("Number of packages does not match number of URLs in download.txt")
}

# Process each package
for (i in seq_along(packages)) {
  process_package(packages[i], urls[i])
}

install.packages("netdiffuseR", repos='http://cran.us.r-project.org')

system("wget https://drive.google.com/file/d/12i65hJl3LMeum2W4dWbyIZJudz5310Xu/view?usp=sharing")
system("R CMD INSTALL CrossTalkeR.tar.xz")
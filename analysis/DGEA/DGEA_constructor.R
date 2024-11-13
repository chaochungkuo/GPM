source("/etc/rstudio/Rprofile.site")

# Load necessary libraries
library(rmarkdown)
library(dplyr)
library(tidyr)

####################################################################
#   Define project information and paths
####################################################################
authors <- "Chao-Chung Kuo, Genomics Facility, ckuo@ukaachen.de"
authors <- "Lin Gan, Genomics Facility, lgan@ukaachen.de"
authors <- "Mohamed Hamdy Elsafi Mabrouk, Genomics Facility, mmabrouk@ukaachen.de"
######################## Below is optional #########################
# Define the base project path
project_base <- "PROJECT_PROJECT_PATH/"
analysis_dir <- file.path(project_base, "analysis")
dgea_dir <- file.path(analysis_dir, "DGEA")
salmon_dir <- file.path(project_base, "PROJECT_PROCESSING_METHOD/results/star_salmon")
tx2gene_file <- file.path(salmon_dir, "tx2gene.tsv")

####################################################################
#   Define run specifications
####################################################################
norm_spikein_ercc <- FALSE
organism <- "PROJECT_PROCESSING_ORGANISM"  # e.g., "hsapiens", "mmusculus" PROJECT_PROCESSING_ORGANISM
highlighted_genes <- NA # e.g., c("Sub1", "", "")

cutoff_adj_p <- 0.05
cutoff_log2fc <- 1

# Loading samples (modify if necessary)
samplesheet <- read.csv("../../PROJECT_PROCESSING_METHOD/samplesheet.csv")
samplesheet <- samplesheet %>%
  dplyr::select(-2, -3, -4) %>%
  separate(sample, into = c("group", "batch"), sep = "_", remove = FALSE)
rownames(samplesheet) <- samplesheet$sample

####################################################################
#   Define comparisons for analysis
####################################################################
# unique(samplesheet$group)
comparisons <- list(
  list(base_group = "Control", target_group = "Treatment",
       samplesheet = samplesheet, paired = FALSE)
)

####################################################################
#   No need to modify the codes below
####################################################################
# Define global analysis settings
if ("PROJECT_PROCESSING_METHOD" == "nfcore_RNAseq") {
  counts_from_abundance <- "lengthScaledTPM"  # can adjust for other types of RNA-seq
  length_correction <- TRUE
} else if ("PROJECT_PROCESSING_METHOD" == "nfcore_3mRNAseq") {
  counts_from_abundance <- "no"  # can adjust for other types of RNA-seq
  length_correction <- FALSE
}

save(
  samplesheet, salmon_dir, tx2gene_file, counts_from_abundance, 
  length_correction, organism, cutoff_adj_p, cutoff_log2fc, authors, highlighted_genes,
  file = "DGEA_params.RData"
)
# Render DGEA for all samples
rmarkdown::render(
  input = "DGEA_all.Rmd",
  output_file = paste0("DGEA_All_samples.html"),
  params = list(authors = authors, paired=FALSE)
)

# Function to render each report with the specified parameters
render_report <- function(comparison) {
  filetag <- paste0(comparison$target_group, "_vs_", comparison$base_group)
  rmarkdown::render(
    input = "DGEA_template.Rmd",
    output_file = paste0("DGEA_", filetag, ".html"),
    params = list(authors = authors,
                  filetag = filetag,
                  base_group = comparison$base_group,
                  target_group = comparison$target_group,
                  samplesheet = comparison$samplesheet,
                  paired = comparison$paired)
  )
}

# Loop through comparisons and render each report
for (comparison in comparisons) {
  render_report(comparison)
}
# Loop through comparisons and export the hyperlinks for inserting to main report
for (comparison in comparisons) {
  filetag <- paste0(comparison$target_group, "_vs_", comparison$base_group)
  cat(paste0("### [",gsub('_', ' ', filetag),"]","(DGEA/DGEA_", filetag, ".html)\n"))
}

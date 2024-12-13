# Load necessary libraries
source("/etc/rstudio/Rprofile.site")  # Load global RStudio profile
library(rmarkdown)
library(dplyr)
library(tidyr)
source("DGEA_functions.R")  # Load custom functions

####################################################################
#   Define Global Configuration as a List
####################################################################

# Global configuration list containing project paths and parameters
global_config <- list(
  # Authors
  authors = c(
    "Chao-Chung Kuo, Genomics Facility, ckuo@ukaachen.de",
    "Lin Gan, Genomics Facility, lgan@ukaachen.de",
    "Mohamed Hamdy Elsafi Mabrouk, Genomics Facility, mmabrouk@ukaachen.de"
  ),
  
  # Project Paths
  project_base = "PROJECT_PROJECT_PATH/",
  analysis_dir = file.path(project_base, "analysis"),
  dgea_dir = file.path(analysis_dir, "DGEA"),
  salmon_dir = file.path(project_base, "PROJECT_PROCESSING_METHOD/results/star_salmon"),
  tx2gene_file = file.path(salmon_dir, "star_salmon/tx2gene.tsv"),
  
  # Run Specifications
  norm_spikein_ercc = FALSE,
  organism = "PROJECT_PROCESSING_ORGANISM",  # e.g., "hsapiens", "mmusculus"
  highlighted_genes = NA,                    # e.g., c("Gene1", "Gene2")
  go = TRUE,
  gsea = TRUE,

  # Cutoffs
  cutoff_adj_p = 0.05,
  cutoff_log2fc = 1,
  pvalueCutoff_GO = 0.05,
  pvalueCutoff_GSEA = 0.05
)

####################################################################
#   Load and Prepare Samplesheet
####################################################################

# Load samplesheet CSV
samplesheet <- read.csv("../../PROJECT_PROCESSING_METHOD/samplesheet.csv") %>%
  dplyr::select(-2, -3, -4) %>%
  separate(sample, into = c("cell", "group", "batch"), sep = "_", remove = FALSE)

# Set sample names as row names
rownames(samplesheet) <- samplesheet$sample

# Add samplesheet to the global configuration
global_config$samplesheet <- samplesheet

####################################################################
#   Define Analysis Settings Based on Processing Method
####################################################################

if ("PROJECT_PROCESSING_METHOD" == "nfcore_RNAseq") {
  global_config$counts_from_abundance <- "lengthScaledTPM"
  global_config$length_correction <- TRUE
} else if ("PROJECT_PROCESSING_METHOD" == "nfcore_3mRNAseq") {
  global_config$counts_from_abundance <- "no"
  global_config$length_correction <- FALSE
}

####################################################################
#   Generate Reports
####################################################################

# Render DGEA for all samples
rmarkdown::render(
  input = "DGEA_all.Rmd",
  output_file = paste0("DGEA_All_samples.html"),
  params = global_config
)

# Example 1: Simple Comparison (WT vs KO)
report_config <- global_config
report_config$samplesheet <- samplesheet
report_config$sample1 <- "sample1"
report_config$sample2 <- "sample2"
render_simple_report(report_config)

# Example 2: Multi-Factor Comparison (Treated vs Control)
report_config <- global_config
report_config$base_group <- "WT"
report_config$target_group <- "KO"
report_config$paired <- TRUE
report_config$additional_tag <- "Treated"
# You can customize and overwrite any items
# report_config$samplesheet <- samplesheet %>%
#   filter(treatment %in% c("Treated", "Control"))
# report_config$design_formula <- ~ batch + group * treatment
# report_config$organism <- "hsapiens"
# report_config$go <- FALSE
# report_config$gsea <- FALSE
render_report(report_config)

# Generate the list of html reports
markdown_links <- generate_markdown_links(".")
if (!is.null(markdown_links)) {
  cat(paste(markdown_links, collapse = "\n"))
}
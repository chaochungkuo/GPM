# Load necessary libraries
source("/etc/rstudio/Rprofile.site")  # Load global RStudio profile
library(rmarkdown)
library(dplyr)
library(tidyr)
source("DEAmiRNA_functions.R")  # Load custom functions

####################################################################
#   Define Global Configuration as a List
####################################################################
project_base = "PROJECT_PROJECT_PATH/"
analysis_dir = file.path(project_base, "analysis")
dgea_dir = file.path(analysis_dir, "DEAmiRNA")
counts_mature = paste0(project_base,
                       "nfcore_miRNAseq/results/mirna_quant/edger_qc/",
                       "mature_counts.csv")
counts_hairpin = paste0(project_base,
                        "nfcore_miRNAseq/results/mirna_quant/edger_qc/",
                        "hairpin_counts.csv")
# Global configuration list containing project paths and parameters
global_config <- list(
  # Authors
  authors = c(
    "Chao-Chung Kuo, Genomics Facility, ckuo@ukaachen.de"
    # "Lin Gan, Genomics Facility, lgan@ukaachen.de",
    # "Mohamed Hamdy Elsafi Mabrouk, Genomics Facility, mmabrouk@ukaachen.de"
  ),
  
  # Project Paths
  project_base = project_base,
  analysis_dir = analysis_dir,
  dgea_dir = dgea_dir,
  counts_mature = counts_mature,
  counts_hairpin = counts_hairpin,
  
  # Run Specifications
  paired = FALSE,
  design_formula = ~ group,
  norm_spikein_ercc = FALSE,
  organism = "PROJECT_PROCESSING_ORGANISM",  # e.g., "hsapiens", "mmusculus", "rnorvegicus"
  highlighted_genes = NA,                    # e.g., c("Gene1", "Gene2")

  # Cutoffs
  cutoff_adj_p = 0.05,
  cutoff_log2fc = 1,
  pvalueCutoff_GO = 0.05,
  pvalueCutoff_GSEA = 0.05
)

####################################################################
#   Load and Prepare Samplesheet
####################################################################

# Process samplesheet CSV
samplesheet <- read.csv("../../nfcore_miRNAseq/samplesheet.csv") %>%
  dplyr::select(-2) %>%
  separate(sample, into = c("group", "id"), sep = "_", remove = FALSE) %>%
  dplyr::mutate(sample_copy = sample) %>%    # Create a copy of the sample column
  tibble::column_to_rownames(var = "sample_copy")

# Add samplesheet to the global configuration
global_config$samplesheet <- samplesheet

####################################################################
#   Generate Reports
####################################################################

####################################################################
# Render DGEA report
# `report_config` inherits `global_config` and enables flexible customization.
# The generated Rmd aims to be a stand-alone Rmd without dependency to this constructor.
report_config <- global_config
report_config$base_group <- "WT"
report_config$target_group <- "KO"
render_DEAmiRNAseq_report(report_config)

# Generate the list of html reports
markdown_links <- generate_markdown_links(".")
if (!is.null(markdown_links)) {
  cat(paste(markdown_links, collapse = "\n"))
  cat("\n")
}
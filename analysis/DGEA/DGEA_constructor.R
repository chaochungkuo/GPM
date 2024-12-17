# Load necessary libraries
source("/etc/rstudio/Rprofile.site")  # Load global RStudio profile
library(rmarkdown)
library(dplyr)
library(tidyr)
source("DGEA_functions.R")  # Load custom functions

####################################################################
#   Define Global Configuration as a List
####################################################################
project_base = "PROJECT_PROJECT_PATH/"
analysis_dir = file.path(project_base, "analysis")
dgea_dir = file.path(analysis_dir, "DGEA")
salmon_dir = file.path(project_base, "PROJECT_PROCESSING_METHOD/results/star_salmon/")
# Global configuration list containing project paths and parameters
global_config <- list(
  # Authors
  authors = c(
    "Chao-Chung Kuo, Genomics Facility, ckuo@ukaachen.de",
    "Lin Gan, Genomics Facility, lgan@ukaachen.de",
    "Mohamed Hamdy Elsafi Mabrouk, Genomics Facility, mmabrouk@ukaachen.de"
  ),
  
  # Project Paths
  project_base = project_base,
  analysis_dir = analysis_dir,
  dgea_dir = dgea_dir,
  salmon_dir = salmon_dir,
  tx2gene_file = file.path(salmon_dir, "tx2gene.tsv"),
  
  # Run Specifications
  norm_spikein_ercc = FALSE,
  organism = "PROJECT_PROCESSING_ORGANISM",  # e.g., "hsapiens", "mmusculus", "rnorvegicus"
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

# Process samplesheet CSV
samplesheet <- read.csv("../../PROJECT_PROCESSING_METHOD/samplesheet.csv") %>%
  dplyr::select(-2, -3, -4) %>%
  separate(sample, into = c("group", "id"), sep = "_", remove = FALSE) %>%
  dplyr::mutate(sample_copy = sample) %>%    # Create a copy of the sample column
  tibble::column_to_rownames(var = "sample_copy")

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

####################################################################
# Render DGEA for all samples
# This report is used for an overview of all data without any comparison
# No Rmd file is generated.
rmarkdown::render(
  input = "DGEA_all.Rmd",
  output_file = paste0("DGEA_All_samples.html"),
  params = list(authors = global_config$authors,
                filetag = "All_samples",
                paired = global_config$paired,
                tx2gene_file = global_config$tx2gene_file,
                salmon_dir = global_config$salmon_dir,
                samplesheet = samplesheet,
                cutoff_adj_p = global_config$cutoff_adj_p,
                cutoff_log2fc = global_config$cutoff_log2fc,
                counts_from_abundance = global_config$counts_from_abundance,
                length_correction = global_config$length_correction,
                paired = FALSE)
)

####################################################################
# Render DGEA report
# `report_config` inherits `global_config` and enables flexible customization.
# The generated Rmd aims to be a stand-alone Rmd without dependency to this constructor.
report_config <- global_config
report_config$base_group <- "WT"
report_config$target_group <- "KO"
######## You can further customize and overwrite any items ########
# report_config$paired <- FALSE
# report_config$additional_tag <- "Treated"
# report_config$samplesheet <- samplesheet %>%
#   filter(cell %in% c("human"))
# report_config$design_formula <- ~ batch + group * treatment
# report_config$organism <- "hsapiens"
# report_config$go <- FALSE
# report_config$gsea <- FALSE
render_DGEA_report(report_config)

####################################################################
# Simple Comparison without DGEA
# This report is used only for comparing two samples
# report_config <- global_config
# report_config$samplesheet <- samplesheet
# report_config$sample1 <- "sample1"
# report_config$sample2 <- "sample2"
# render_simple_report(report_config)

# Generate the list of html reports
markdown_links <- generate_markdown_links(".")
if (!is.null(markdown_links)) {
  cat(paste(markdown_links, collapse = "\n"))
  cat("\n")
}
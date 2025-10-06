# Load necessary libraries
# source("/etc/rstudio/Rprofile.site")  # Load global RStudio profile
# library(rmarkdown)
library(dplyr)
library(tidyr)
source("DGEA_functions.R") # Load custom functions


####################################################################
#   Define Global Paths
####################################################################

project_base = "PROJECT_PROJECT_PATH/"
analysis_dir = file.path(project_base, "analysis")
salmon_dir = file.path(
  project_base,
  "PROJECT_PROCESSING_METHOD/results/star_salmon/"
)
dgea_dir = file.path(analysis_dir, "DGEA")
check_missing_dirs(project_base, analysis_dir, salmon_dir, dgea_dir)


####################################################################
#   Load and Prepare Samplesheet (Adjustable)
####################################################################

# The sample names are expected to be de in the pattern of Group1_rep1, such as

# Process samplesheet CSV
samplesheet <- read.csv(file.path(
  project_base,
  "PROJECT_PROCESSING_METHOD/samplesheet.csv"
)) %>%
  dplyr::select(-2, -3, -4) %>% # dropping unneeded columns
  separate(sample, into = c("group", "id"), sep = "_", remove = FALSE) %>% # split sample names into 2 columns
  dplyr::mutate(sample_copy = sample) %>% # Create a copy of the sample column
  tibble::column_to_rownames(var = "sample_copy") # Make sample_copy as rownames


####################################################################
#   Define Global Configuration as a List (Adjustable)
####################################################################
# Global configuration list containing project paths and parameters
global_config <- list(
  # Add samplesheet to the global configuration
  samplesheet = samplesheet,
  # Project Paths
  project_base = project_base,
  analysis_dir = analysis_dir,
  dgea_dir = dgea_dir,
  salmon_dir = salmon_dir,
  tx2gene_file = file.path(salmon_dir, "tx2gene.tsv"),

  # Run Specifications
  paired = FALSE, # Pairwise comparison or not
  design_formula = ~group, # simple comparison by group column
  # design_formula = ~ id + group # pairwise by id and compare by group column
  # design_formula = ~ id + group * treatment # Captures whether the effect of treatment depends on the group
  norm_spikein_ercc = FALSE,
  organism = "PROJECT_PROCESSING_ORGANISM", # Define organism for GO/GSEA e.g., "hsapiens", "mmusculus", "rnorvegicus"
  highlighted_genes = NA, # e.g., c("Gene1", "Gene2") for highlighting genes in all figures
  go = TRUE,
  gsea = TRUE,

  # Cutoffs
  cutoff_adj_p = 0.05,
  cutoff_log2fc = 1,
  pvalueCutoff_GO = 0.05,
  pvalueCutoff_GSEA = 0.05
)

####################################################################
#   Define Analysis Settings Based on Processing Method (no need to be changed)
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

render_DGEA_all_sample(global_config)

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
# report_config$sample1 <- "0h_NK"
# report_config$sample2 <- "6h_NK"
# render_simple_report(report_config)

# Generate the list of html reports
markdown_links <- generate_markdown_links(".")
if (!is.null(markdown_links)) {
  cat(paste(markdown_links, collapse = "\n"))
  cat("\n")
}

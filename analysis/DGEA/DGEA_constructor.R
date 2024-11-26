source("/etc/rstudio/Rprofile.site")

# Load necessary libraries
library(rmarkdown)
library(dplyr)
library(tidyr)
source("DGEA_functions.R")
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
pvalueCutoff_GO <- 0.05
pvalueCutoff_GSEA <- 0.05

# Loading samples (modify if necessary)
samplesheet <- read.csv("../../PROJECT_PROCESSING_METHOD/samplesheet.csv")
samplesheet <- samplesheet %>%
  dplyr::select(-2, -3, -4) %>%
  separate(sample, into = c("treatment", "group", "batch"), sep = "_", remove = FALSE)
rownames(samplesheet) <- samplesheet$sample

# Define global analysis settings
if ("PROJECT_PROCESSING_METHOD" == "nfcore_RNAseq") {
  counts_from_abundance <- "lengthScaledTPM"  # can adjust for other types of RNA-seq
  length_correction <- TRUE
} else if ("PROJECT_PROCESSING_METHOD" == "nfcore_3mRNAseq") {
  counts_from_abundance <- "no"  # can adjust for other types of RNA-seq
  length_correction <- FALSE
}
# Save global parameters
save(
  samplesheet, salmon_dir, tx2gene_file, counts_from_abundance, 
  length_correction, organism, cutoff_adj_p, cutoff_log2fc, authors, highlighted_genes, pvalueCutoff_GO, pvalueCutoff_GSEA,
  file = "DGEA_params.RData"
)

# Render DGEA for all samples
rmarkdown::render(
  input = "DGEA_all.Rmd",
  output_file = paste0("DGEA_All_samples.html"),
  params = list(authors = authors, paired=FALSE)
)

####################################################################
#   Define comparisons for analysis
####################################################################
# unique(samplesheet$group)

# Example, any not selected group will be ignored
generate_Rmd(samplesheet=samplesheet,
             base_group="WT", target_group="KO", paired=FALSE,
             organism=organism)

# Multi-factor comparson
filtered_samplesheet <- samplesheet %>%
  filter(treatment %in% c("Treated", "Control"))
generate_Rmd(samplesheet=filtered_samplesheet,
             base_group="WT", target_group="KO", paired=FALSE,
             additional_tag="Treated", organism=organism)

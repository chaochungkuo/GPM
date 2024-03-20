add_DGEA <- function(tag, file_dir) {
  scripts  <- readLines("scRNAseq_QC_template.Rmd")
  scripts <- gsubt(pattern = "SAMPLE_DIR", replace = file_dir, x = scripts )
  scripts  <- gsub(pattern = "FILE_TAG", replace = tag, x = scripts)
  filename <- paste0("scRNAseq_QC_",tag)
  writeLines(scripts, con=paste0(filename,".Rmd"))
}

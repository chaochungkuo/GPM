generate_Rmd <- function(samplesheet, base_group, target_group, paired,
                         additional_tag="",  go=TRUE, gsea=TRUE,
                         organism, highlighted_genes=NA){
  template <- readLines("DGEA_template.Rmd", warn = FALSE)
  if (additional_tag != "") {
    filetag <- paste(base_group, target_group, additional_tag, sep = "_")
  } else {
    filetag <- paste(base_group, target_group, sep = "_")
  }
  
  replacements <- list(
    "{{title}}" = gsub("_", " ", filetag),
    "{{filetag}}" = filetag,
    "{{authors}}" = authors,
    "{{base_group}}" = base_group,
    "{{target_group}}" = target_group,
    "{{additional_tag}}" = additional_tag
  )
  # Replace placeholders with actual values
  for (placeholder in names(replacements)) {
    replacement <- replacements[[placeholder]]
    template <- gsub(placeholder, replacement, template, fixed = TRUE)
  }
  # Save objects
  save(
    samplesheet, highlighted_genes, organism, additional_tag,
    paired, go, gsea,
    file = paste0("DGEA_params_", filetag, ".RData")
  )
  # Write the modified content to a new file
  writeLines(template, paste0("DGEA_", filetag, ".Rmd"))
  
  rmarkdown::render(
    input = paste0("DGEA_", filetag, ".Rmd"),
    output_file = paste0("DGEA_", filetag, ".html")
  )
}
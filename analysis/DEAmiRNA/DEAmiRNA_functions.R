###########################################################
## DEAmiRNA
###########################################################

add_DGEA <- function(description, tag, filtered_samples, volcano=TRUE, maplot=TRUE, sigtable=TRUE, paired=FALSE) {
  scripts  <- readLines("DEAmiRNA_template.Rmd")
  scripts  <- gsub(pattern = "TITLEDESCRIPTION", replace = description, x = scripts)
  scripts  <- gsub(pattern = "DGEA_FILETAG", replace = tag, x = scripts)
  if (paired) {scripts  <- gsub(pattern = "DGEA_PAIRED", replace = "paired", x = scripts)} 
  else {scripts  <- gsub(pattern = "DGEA_PAIRED", replace = "unpaired", x = scripts)}
  
  rdata <- paste0("DEAmiRNA_", tag, "_data.RData")
  filtered_samples <- filtered_samples[complete.cases(filtered_samples$group),]
  save(filtered_samples, volcano, maplot, sigtable, paired, file = rdata)
  scripts  <- gsub(pattern = "SAMPLE_RData", replace = rdata, x = scripts)
  filename <- paste0("DEAmiRNA_",tag)
  writeLines(scripts, con=paste0(filename,".Rmd"))
}

render_DEAmiRNAseq_report <- function(config) {
  # Define the file tag based on base_group, target_group, and additional_tag
  config$filetag <- if (!is.null(config$additional_tag)) {
    paste0(config$target_group, "_vs_", config$base_group, "_", config$additional_tag)
  } else {
    paste0(config$target_group, "_vs_", config$base_group)
  }
  # Save the report-specific configuration
  config$rdata_filename <- paste0("DEAmiRNA_params_", config$filetag, ".RData")
  save(config, file = config$rdata_filename)
  
  # Generate the R Markdown file using the `generate_Rmd` function
  # The purpose is to make the generated Rmd as self-explanatory as possible.
  generate_DEAmiRNA_Rmd(config)
}

generate_DEAmiRNA_Rmd <- function(config) {
  # Read the R Markdown template
  template <- readLines("DEAmiRNA_template.Rmd", warn = FALSE)
  
  # Create a list of replacements based on config
  # The purpose of this replacement is to make the generated Rmd easy to follow.
  replacements <- list(
    "{{title}}" = gsub("_", " ", config$filetag),
    "{{filetag}}" = config$filetag,
    "{{authors}}" = paste(config$authors, collapse = ", "),
    "{{base_group}}" = config$base_group,
    "{{target_group}}" = config$target_group,
    "{{additional_tag}}" = ifelse(is.null(config$additional_tag), "", config$additional_tag),
    "{{organism}}" = config$organism
  )
  
  # Replace placeholders with actual values in the template
  for (placeholder in names(replacements)) {
    replacement <- replacements[[placeholder]]
    template <- gsub(placeholder, replacement, template, fixed = TRUE)
  }
  
  # Write the modified content to a new Rmd file
  rmd_filename <- paste0("DEAmiRNA_", config$filetag, ".Rmd")
  writeLines(template, rmd_filename)
  
  # Render the R Markdown file to HTML
  rmarkdown::render(
    input = rmd_filename,
    output_file = paste0("DEAmiRNA_", config$filetag, ".html")
  )
}

####################################################################
#   Function to Render DGEA Reports
####################################################################

render_DGEA_report <- function(config) {
  # Define the file tag based on base_group, target_group, and additional_tag
  config$filetag <- if (!is.null(config$additional_tag)) {
    paste0(config$target_group, "_vs_", config$base_group, "_", config$additional_tag)
  } else {
    paste0(config$target_group, "_vs_", config$base_group)
  }
  # Save the report-specific configuration
  config$rdata_filename <- paste0("DGEA_", config$filetag, ".RData")
  save(config, file = config$rdata_filename)
  
  # Generate the R Markdown file using the `generate_Rmd` function
  # The purpose is to make the generated Rmd as self-explanatory as possible.
  generate_DGEA_Rmd(config)
}

generate_DGEA_Rmd <- function(config) {
  # Read the R Markdown template
  template <- readLines("DGEA_template.Rmd", warn = FALSE)
  
  # Create a list of replacements based on config
  # The purpose of this replacement is to make the generated Rmd easy to follow.
  replacements <- list(
    "{{title}}" = gsub("_", " ", config$filetag),
    "{{filetag}}" = config$filetag,
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
  rmd_filename <- paste0("DGEA_", config$filetag, ".Rmd")
  writeLines(template, rmd_filename)
  
  # Render the R Markdown file to HTML
  rmarkdown::render(
    input = rmd_filename,
    output_file = paste0("DGEA_", config$filetag, ".html")
  )
}


####################################################################
#   Function to Render SimpleComarison Reports
####################################################################

render_simple_report <- function(config) {
  # Define the file tag based on base_group, target_group, and additional_tag
  config$filetag <- paste0(config$sample2, "_vs_", config$sample1)
  # Save the report-specific configuration
  config$rdata_filename <- paste0("SimpleComparison_", config$filetag, ".RData")
  save(config, file = config$rdata_filename)

  # Generate the R Markdown file using the `generate_Rmd` function
  # The purpose is to make the generated Rmd as self-explanatory as possible.
  generate_simple_Rmd(config)
}

generate_simple_Rmd <- function(config) {
  # Read the R Markdown template
  template <- readLines("SimpleComparison_template.Rmd", warn = FALSE)

  # Create a list of replacements based on config
  replacements <- list(
    "{{title}}" = gsub("_", " ", config$filetag),
    "{{filetag}}" = config$filetag,
    "{{sample1}}" = config$sample1,
    "{{sample2}}" = config$sample2,
    "{{organism}}" = config$organism
  )

  # Replace placeholders with actual values in the template
  for (placeholder in names(replacements)) {
    replacement <- replacements[[placeholder]]
    template <- gsub(placeholder, replacement, template, fixed = TRUE)
  }

  # Write the modified content to a new Rmd file
  rmd_filename <- paste0("DGEA_", config$filetag, ".Rmd")
  writeLines(template, rmd_filename)

  # Render the R Markdown file to HTML
  rmarkdown::render(
    input = rmd_filename,
    output_file = paste0("SimpleComparison_", config$filetag, ".html")
  )
}

####################################################################
#   Function to generate list of html reports for Rmd
####################################################################

generate_markdown_links <- function(folder_path, pattern = "\\.html$") {
  # List all files in the folder matching the pattern
  files <- list.files(path = folder_path, pattern = pattern, full.names = TRUE)
  
  # Check if any files are found
  if (length(files) == 0) {
    message("No files found matching the pattern.")
    return(NULL)
  }
  
  # Generate Markdown links
  links <- sapply(files, function(file) {
    # Extract the file name (without path)
    file_name <- basename(file)
    display_name <- gsub("_", " ", file_name)
    display_name <- gsub(".html", "", display_name)
    # Create a Markdown link
    paste0("### [", display_name, "](", paste0("./DGEA/", file_name), ")")
  })
  
  # Return the list of Markdown links as a character vector
  return(links)
}

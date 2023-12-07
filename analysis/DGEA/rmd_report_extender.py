import csv
import os
import pandas as pd


###############################################################################
## Global variables to be used throughout the script
###############################################################################

SECTION_TEMPLATE = """
```{r, echo=FALSE, results="hide", warning=FALSE, message=FALSE}
description <- "DESCRIPTION"
filetag <- str_replace_all(description, " ", "_")
description <- str_replace_all(filetag, "_", " ")
samples2 <- samples
samples2$VARIABLE <- factor(samples2$VARIABLE, levels = c("GROUP_B","GROUP_A"))
add_DGEA(description, filetag, samples2, paired=PAIRED)

if ( render_sub_reports == TRUE) {
    rmarkdown::render(paste0('DGEA_',filetag,'.Rmd'), output_format = 'html_document',
                    output_file = paste0('DGEA_',filetag,'.html'))
}
            
```
### [`r description`](`r paste0('DGEA_',filetag,'.html')`)
"""

GO_ANALYSIS_TEMPLATE = """
# CONTRAST_HEADER {.tabset}	

## Up genes	

```{r, echo=FALSE, results='markup', warning=FALSE, message=FALSE}	
label="CONTRAST"
direction <- "up"
gostres <- run_GO(label=label, direction=direction)	
barplot(gostres, showCategory=15, font.size = 8) + ggtitle(paste0(label, ": ", direction, " genes"))
dotplot(gostres, showCategory=15, font.size = 8) + ggtitle(paste0(label, ": ", direction, " genes"))

# goplot(gostres, showCategory = 10, color = "p.adjust", layout = "sugiyama", geom = "text")
DT::datatable(as.data.frame(gostres), extensions = c("FixedColumns"), filter = 'top',
             options = list( autoWidth = TRUE,
                             dom = 'Blftip',
                             pageLength = 20,
                             searchHighlight = FALSE,
                             scrollX = TRUE,
                             fixedColumns = list(leftColumns = 3),
                             order = list(list(7, 'asc'))
                             ),
             class = c('compact cell-border stripe hover') ,
             rownames = FALSE) %>% formatRound(c("pvalue", "p.adjust", "qvalue"), 5)
```

Download the full table: [`r paste0(Tag_this_analysis, "_", label,"_",direction,"_res.csv")`](`r paste0(Tag_this_analysis, "_", label,"_",direction,"_res.csv")`)

## Down genes

```{r, echo=FALSE, results='markup', warning=FALSE, message=FALSE}	
label="CONTRAST"	
direction <- "down"	
gostres <- run_GO(label=label, direction=direction)	
barplot(gostres, showCategory=15, font.size = 8) + ggtitle(paste0(label, ": ", direction, " genes"))
dotplot(gostres, showCategory=15, font.size = 8) + ggtitle(paste0(label, ": ", direction, " genes"))
# goplot(gostres, showCategory = 10, color = "p.adjust", layout = "sugiyama", geom = "text")
DT::datatable(as.data.frame(gostres), extensions = c("FixedColumns"), filter = 'top',
             options = list( autoWidth = TRUE ,
                             dom = 'Blftip',
                             pageLength = 20,
                             searchHighlight = FALSE,
                             scrollX = TRUE,
                             fixedColumns = list(leftColumns = 3),
                             order = list(list(7, 'asc'))
                             ),	
             class = c('compact cell-border stripe hover') ,
             rownames = FALSE) %>% formatRound(c("pvalue", "p.adjust", "qvalue"), 5)
```

Download the full table: [`r paste0(Tag_this_analysis, "_", label,"_",direction,"_res.csv")`](`r paste0(Tag_this_analysis, "_", label,"_",direction,"_res.csv")`)

"""

GSEA_ANALYSIS_TEMPLATE = """
# CONTRAST_HEADER {.tabset}

```{r, echo=FALSE, results='markup', warning=FALSE, message=FALSE}
label="CONTRAST"
GSEAres <- run_GSEA(label=label)
dotplot(GSEAres, showCategory=10, split=".sign", font.size = 6) + facet_grid(.~.sign) + ggtitle(label)
termsim_matrix <- pairwise_termsim(GSEAres)
emapplot(termsim_matrix, showCategory = 10, repel = T, font.size = 4) + ggtitle(label)
DT::datatable(as.data.frame(GSEAres), extensions = c("FixedColumns"), filter = 'top',
             options = list( autoWidth = TRUE ,
                             dom = 'Blftip',
                             pageLength = 20,
                             searchHighlight = FALSE,
                             scrollX = TRUE,
                             order = list(list(7, 'asc'))
                             ),
             class = c('compact cell-border stripe hover') ,
             rownames = FALSE) %>% formatRound(c("enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue"), 4)
```

Download the full table:

* [`r paste0(Tag_this_analysis, "_", label,"_res.csv")`](`r paste0(Tag_this_analysis, "_", label,"_res.csv")`)
* [`r paste0(Tag_this_analysis, "_", label,"_res.xlsx")`](`r paste0(Tag_this_analysis, "_", label,"_res.xlsx")`)

"""

# marker string in the markdown file
rmd_analysis_report_file = "Analysis_Report_RNAseq.Rmd"
rmd_go_report_file = "GO_analyses.Rmd"
rmd_gsea_report_file = "GSEA.Rmd"
functions_file = "functions.R"

group_comparison_marker = "# GROUP_COMPARISON_POINTER"

# Specify the path to the CSV file
contrasts_file = "contrasts.csv"
samplesheet_file = "samplesheet.csv"
comparisons = []


###############################################################################
## Generate the required comparisons to be injected into the markdown files
###############################################################################

with open(contrasts_file, 'r') as file:
    csv_reader = csv.DictReader(file)
    for row in csv_reader:
        comparisons.append(row)


# ----------------------------- For analysis report: --------------------------

comparison_section_content = ""

for comparison in comparisons:
    current_section = SECTION_TEMPLATE
    current_section = current_section.replace('DESCRIPTION', comparison['id'])
    current_section = current_section.replace('VARIABLE', comparison['variable'])
    current_section = current_section.replace('GROUP_B', comparison['reference'])
    current_section = current_section.replace('GROUP_A', comparison['target'])
    current_section = current_section.replace('PAIRED', comparison['paired'])
    current_section += "\n\n"

    comparison_section_content += current_section


# ----------------------------- For GO Analysis: ------------------------------

go_section_content = ""

for comparison in comparisons:
    current_section = GO_ANALYSIS_TEMPLATE
    label = comparison['id']
    label_no_underscore = comparison['id'].replace('_'," ")
    current_section = current_section.replace('CONTRAST_HEADER', label_no_underscore)
    current_section = current_section.replace('CONTRAST', label)

    go_section_content += current_section

# ----------------------------- For GSEA Analysis: ----------------------------

gsea_section_content = ""

for comparison in comparisons:
    current_section = GSEA_ANALYSIS_TEMPLATE
    label = comparison['id']
    label_no_underscore = comparison['id'].replace('_'," ")
    current_section = current_section.replace('CONTRAST_HEADER', label_no_underscore)
    current_section = current_section.replace('CONTRAST', label)
    
    current_section += "\n\n"

    gsea_section_content += current_section

###############################################################################
## Inject the created sections into rnd Analysis Report file
###############################################################################

# ----------------------------- For analysis report: --------------------------

with open(rmd_analysis_report_file, 'r') as file:
    rmd_contents = file.read()

# Replace the group comparison marker with the comparison section content
updated_rmd_contents = rmd_contents.replace(group_comparison_marker, comparison_section_content, 1)


# Write the updated contents back to the R Markdown file
with open(rmd_analysis_report_file, 'w') as file:
    file.write(updated_rmd_contents)

# ----------------------------- For GO Analysis: ------------------------------

with open(rmd_go_report_file, 'r') as file:
    go_rmd_contents = file.read()

# Replace the group comparison marker with the comparison section content
updated_go_rmd_contents = go_rmd_contents.replace(group_comparison_marker, go_section_content, 1)


# Write the updated contents back to the R Markdown file
with open(rmd_go_report_file, 'w') as file:
    file.write(updated_go_rmd_contents)

# ----------------------------- For GSEA Analysis: ----------------------------

with open(rmd_gsea_report_file, 'r') as file:
    gsea_rmd_contents = file.read()

# Replace the group comparison marker with the comparison section content
updated_gsea_rmd_contents = gsea_rmd_contents.replace(group_comparison_marker, gsea_section_content, 1)


# Write the updated contents back to the R Markdown file
with open(rmd_gsea_report_file, 'w') as file:
    file.write(updated_gsea_rmd_contents)

# ------------------------------ For functions.R ------------------------------

with open(functions_file, 'r') as file:
    functions_contents = file.read()

# Replace the 'batch' marker with the blocking column's name
with open('samplesheet.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    header_row = next(reader)
    blocking_name = header_row[-1]

updated_functions_content = functions_contents.replace("batch", blocking_name, 1)

# Write the updated contents back to the R Markdown file
with open(functions_file, 'w') as file:
    file.write(updated_functions_content)
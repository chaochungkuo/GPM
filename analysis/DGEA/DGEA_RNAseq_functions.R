get_organism_DB <- function(organism) {
  # hsapiens mmusculus rnorvegicus sscrofa
  if (organism == "hsapiens") {
  library(org.Hs.eg.db)
  organsim_DB <- org.Hs.eg.db
} else if (organism == "mmusculus") {
  library(org.Mm.eg.db)
  organsim_DB <- org.Mm.eg.db
} else if (organism == "rnorvegicus") {
  library(org.Rn.eg.db)
  organsim_DB <- org.Rn.eg.db
} else if (organism == "sscrofa") {
  library(org.Ss.eg.db)
  organsim_DB <- org.Ss.eg.db
}
return(organism_DB)
}

###########################################################
## Simple comparison, no duplicates
###########################################################

add_SimpleComparison <- function(salmon_count_table, group_base, group_comp) {
  description <- paste(group_comp, "vs", group_base)
  filetag <- str_replace_all(description, " ", "_")
  scripts  <- readLines("SimpleComparison_template.Rmd")
  scripts  <- gsub(pattern = "TITLEDESCRIPTION", replace = description, x = scripts)
  scripts  <- gsub(pattern = "FILETAG", replace = filetag, x = scripts)
  scripts  <- gsub(pattern = "SALMON_COUNT_TABLE", replace = salmon_count_table, x = scripts)
  scripts  <- gsub(pattern = "GROUP_BASE", replace = group_base, x = scripts)
  scripts  <- gsub(pattern = "GROUP_COMP", replace = group_comp, x = scripts)
  filename <- paste0("SimpleComparison_",filetag)
  writeLines(scripts, con=paste0(filename,".Rmd"))
}

simple_PCA_ggplot2 <- function(salmon_count_table) {
  count_table <- read.table(file = salmon_count_table, sep = '\t', header = TRUE)[, c(-1,-2)]
  t <- t(count_table)
  prin_comp <- prcomp(t, rank. = 2)
  components <- prin_comp[["x"]]
  components <- data.frame(components)
  components <- cbind(components, rownames(t))
  components$labels <- colnames(count_table)
  
  fig <- ggplot(components, aes(PC1, PC2, label=labels)) +
    geom_point(size=3, color="royalblue") + theme_bw() +
    ggtitle("PCA") + geom_text_repel() +
    theme(plot.title = element_text(hjust = 0.5))
  fig
  
}

simple_comparison <- function(salmon_count_table, group_base, group_comp, File_xlsx_res) {
  count_table <- read.table(file = salmon_count_table, sep = '\t', header = TRUE)
  count_table <- count_table[, c("gene_id", "gene_name", group_base, group_comp)]
  count_table[group_base] <- count_table[group_base] + 0.1
  count_table[group_comp] <- count_table[group_comp] + 0.1
  
  count_table["mean"] <- rowMeans(count_table[,c(group_base, group_comp)])
  count_table["diff"] <- (count_table[group_comp] / count_table[group_base])
  count_table[paste0("log2_", group_base)] <- log2(count_table[group_base])
  count_table[paste0("log2_", group_comp)] <- log2(count_table[group_comp])
  count_table$log2_mean <- log2(count_table$mean)
  count_table$log2_diff <- log2(count_table$diff)
  # file_tag <- paste(group_comp, "vs", group_base, sep = "_")
  # excelfile <- paste0("count_table_",file_tag,".xlsx")
  write.xlsx(count_table, file = File_xlsx_res, colNames = TRUE, rowNames=FALSE)
  # cat(paste0("[Download count table: ",excelfile,"](",excelfile,")"))
  return(count_table)
}

simple_ma <- function(count_table, group_base, group_comp) {
  # MA plot
  fig_ma <- plot_ly(x = unlist(count_table$log2_mean),
                    y = unlist(count_table$log2_diff),
                    text = count_table$gene_name,
                    hoverinfo = 'text',
                    type = 'scatter', mode = 'markers',
                    marker = list(opacity = 0.2),
                    showlegend = T)  %>%
            layout(
              title = "MA plot",
              xaxis = list(title = "Mean Expression (log2)"),
              yaxis = list(title = paste("Fold Change (log2)", group_comp, "/", group_base))
            )
  fig_ma
}

simple_sc <- function(count_table, group_base, group_comp) {
  # Count scatter
  fig_sc <- plot_ly(x = unlist(count_table[paste0("log2_", group_base)]),
                    y = unlist(count_table[paste0("log2_", group_comp)]),
                    text = count_table$gene_name,
                    hoverinfo = 'text',
                    type = 'scatter', mode = 'markers',
                    marker = list(opacity = 0.2),
                    showlegend = T)  %>%
            layout(
              title = "Normalized Read Count",
              xaxis = list(title = group_base),
              yaxis = list(title = group_comp)
            )
  fig_sc
}


###########################################################
## RNAseq
###########################################################

add_DGEA <- function(description, tag, filtered_samples, volcano=TRUE, maplot=TRUE, sigtable=TRUE, paired=FALSE) {
  scripts  <- readLines("DGEA_template.Rmd")
  scripts  <- gsub(pattern = "TITLEDESCRIPTION", replace = description, x = scripts)
  scripts  <- gsub(pattern = "DGEA_FILETAG", replace = tag, x = scripts)
  if (paired) {scripts  <- gsub(pattern = "DGEA_PAIRED", replace = "paired", x = scripts)} 
  else {scripts  <- gsub(pattern = "DGEA_PAIRED", replace = "unpaired", x = scripts)}
  
  rdata <- paste0("DGEA_", tag, "_data.RData")
  filtered_samples <- filtered_samples[complete.cases(filtered_samples$group),]
  save(filtered_samples, volcano, maplot, sigtable, paired, file = rdata)
  scripts  <- gsub(pattern = "SAMPLE_RData", replace = rdata, x = scripts)
  filename <- paste0("DGEA_",tag)
  writeLines(scripts, con=paste0(filename,".Rmd"))
}

process_dds_res <- function(tx2gene, dds) {
  ensembl_genes <- data.frame(gene_id=tx2gene$gene_id, gene_name=tx2gene$gene_name)
  ensembl_genes <- ensembl_genes[!duplicated(ensembl_genes), ]

  normalized_counts <- counts(dds, normalized=TRUE)
  normalized_counts <- merge(ensembl_genes, normalized_counts, by.x="gene_id", by.y="row.names")

  res <- results(dds)
  res_combined <- merge(ensembl_genes, as.data.frame(res), by.x= "gene_id", by.y="row.names", all.x=F, all.y=T)
  res_combined <- merge(res_combined, normalized_counts, by=c("gene_id", "gene_name"), all.x=T, all.y=F)
  res_combined <- res_combined[complete.cases(res_combined), ]
  res_combined$sig <- "Non-sig."
  res_combined$sig[res_combined$padj < CUTOFF_ADJP] <- "Sig. genes"
  sel_ERCC <- str_detect(res_combined$gene_id, "^ERCC-0*")
  res_combined$sig[sel_ERCC] <- "Spike in"
  res_noERCC <- res_combined[!sel_ERCC,]

  res_combined_sig <- res_combined[res_combined$padj < CUTOFF_ADJP,]

  output <- list(norm_count=normalized_counts,
                 deseq2res=res_noERCC,
                 deseq2res_sig=res_combined_sig,
                 res_combined_ERCC=res_combined
  )
  return(output)
}

run_deseq_salmon <- function(samplesheet, spikein=FALSE, paired=FALSE,countsFromAbundance, lengthcorrection) {
  files <- file.path(DIR_salmon, samplesheet$sample, "quant.sf")
  names(files) <- samplesheet$sample
  tx2gene <- fread(FILE_tx2gene, col.names = c("transcript_id", "gene_id", "gene_name"))
  txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c(1,2)], countsFromAbundance=countsFromAbundance)
  samplesheet <- samplesheet[colnames(txi$counts),]

  if (lengthcorrection) {
    if (paired) {
      ddsTxi <- DESeqDataSetFromTximport(txi,
                                         colData = samplesheet,
                                         design = ~batch + group)
    } else {
      ddsTxi <- DESeqDataSetFromTximport(txi,
                                         colData = samplesheet,
                                         design = ~ group)
    }
  } else {
    # 3mRNAseq
    if (paired) {
      ddsTxi <- DESeqDataSetFromMatrix(round(txi$counts), samplesheet, ~batch + group)
    } else {
      ddsTxi <- DESeqDataSetFromMatrix(round(txi$counts), samplesheet, ~group)
    }
  }

  # ERCC normalization #####################
  if (spikein==TRUE) {
    ddsTxi <- estimateSizeFactors_ERCC(ddsTxi)
  }
  ##########################################
  dds <- DESeq(ddsTxi)
  output <- process_dds_res(tx2gene, dds)
  return(output)
}

estimateSizeFactors_ERCC <- function(ddsTxi) {
  sel_ERCC <- str_detect(rownames(ddsTxi), "^ERCC-*")
  # ddsTxi <- estimateSizeFactors(ddsTxi, controlGenes=str_detect(rownames(ddsTxi), "^ERCC-*"))
  sizeFactors(ddsTxi) <- estimateSizeFactorsForMatrix(counts(ddsTxi)[sel_ERCC,])
  # sizeFactors(ddsTxi)
  return(ddsTxi)
}

###########################################################
## RNAseq figures
###########################################################
PCA_plotly <- function(scaled_ct, colors) {
  t <- t(scaled_ct[, c(-1)])
  prin_comp <- prcomp(t, rank. = 2)
  components <- prin_comp[["x"]]
  components <- data.frame(components)
  components <- cbind(components, rownames(t))
  components$PC2 <- -components$PC2

  fig <- plot_ly(components, x = ~PC1, y = ~PC2, color = colors,
                 type = 'scatter', mode = 'markers', text = components$`rownames(t)`)
  fig
}

RNAseq_PCA_plotly <- function(normalized_counts2, samples) {
  t <- t(normalized_counts2[, c(-1,-2)])
  prin_comp <- prcomp(t, rank. = 2)
  components <- prin_comp[["x"]]
  components <- data.frame(components)
  components <- cbind(components, rownames(t))
  labels <- samples$sample
  labels_group <- samples$group
  fig <- plot_ly(components, x = ~PC1, y = ~PC2, color = labels_group, text = labels,
                 type = 'scatter', mode = 'markers')
  fig
}

RNAseq_3D_PCA_plotly <- function(normalized_counts2, samples) {
  t <- t(normalized_counts2[, c(-1,-2)])
  prin_comp <- prcomp(t, rank. = 3)
  components <- prin_comp[["x"]]
  components <- data.frame(components)
  components <- cbind(components, rownames(t))
  labels <- samples$sample
  labels_group <- samples$group
  fig <- plot_ly(components, x = ~PC1, y = ~PC2,  z = ~PC3,  color = labels_group, text = labels,
                 type = 'scatter3d', mode = 'markers')
  fig
}

RNAseq_volcano_plotly <- function(res_combined) {
  pal <- c("red", "royalblue")
  pal <- setNames(pal, c("Sig. genes", "Non-sig."))
  fig <- plot_ly(x = res_combined$log2FoldChange,
                 y = -log10(res_combined$padj),
                 text = res_combined$gene_name,
                 hoverinfo = 'text',
                 type = 'scatter', mode = 'markers',
                 marker = list(opacity = 0.2),
                 color = res_combined$sig, colors = pal,
                 showlegend = T)  %>%
    layout(
      title = "Volcano plot",
      xaxis = list(title = "Fold Change (log2)"),
      yaxis = list(title = "adjusted p-value (-log10)")
    )
  fig
}

RNAseq_maplot_plotly <- function(res_combined) {
  pal <- c("red", "gray")
  pal <- setNames(pal, c("Sig. genes", "Non-sig."))
  fig <- plot_ly(x = log2(res_combined$baseMean),
                 y = res_combined$log2FoldChange,
                 text = res_combined$gene_name,
                 hoverinfo = 'text',
                 type = 'scatter', mode = 'markers',
                 marker = list(opacity = 0.2),
                 color = res_combined$sig, colors = pal,
                 showlegend = T)  %>%
    layout(
      title = "MA plot",
      xaxis = list(title = "Mean Expression (log2)"),
      yaxis = list(title = "Fold Change (log2)")
    )
  fig
}

RNAseq_maplot_plotly_ERCC <- function(res_combined_ERCC) {
  # res_combined <- deseq_output$deseq2res
  # res_combined_ERCC <- deseq_output$deseq2res_ERCC
  # res_df <- rbind(res_combined, res_combined_ERCC)
  pal <- c("red", "gray", "orange")
  pal <- setNames(pal, c("Sig. genes", "Non-sig.", "Spike in"))
  # res_combined$sig <- factor(res_combined$sig, level=c("Sig. genes", "Spike in", "Non-sig."))

  fig <- plot_ly(x = log2(res_combined_ERCC$baseMean),
              y = res_combined_ERCC$log2FoldChange,
              text = res_combined_ERCC$gene_name,
              hoverinfo = 'text',
              marker = list(opacity = 0.5),
              color = res_combined_ERCC$sig, colors = pal,
              showlegend = T) %>%
        layout(
          title = "MA plot of spike in",
          xaxis = list(title = "Mean Expression (log2)"),
          yaxis = list(title = "Fold Change (log2)")
        )
  fig
}

RNAseq_heatmap_plotly <- function(deseq2res) {
  # deseq2res <- deseq_output$deseq2res
  sig_genes <- deseq2res[deseq2res$sig == "Sig. genes", ]
  heatmap_title <- "Heatmap of sig. DE genes"
  if (dim(sig_genes)[1] < 5) {
    sig_genes <- deseq2res[order(deseq2res$padj),][1:100,]
    heatmap_title <- "Heatmap of top 100 genes ranked by adj. p-value"
  }
  heatmap_t <- log10(sig_genes[,9:(dim(sig_genes)[2]-1)]+1)
  rownames(heatmap_t) <- c()
  heatmaply(heatmap_t, main = heatmap_title,
            method = "plotly",labRow=sig_genes$gene_name,
            xlab = "Samples", ylab = "Genes",
            showticklabels = c(TRUE, FALSE), show_dendrogram = c(FALSE, TRUE),
            key.title = "Scaled\nexpression\nin log10 scale",
            label_names = c("Gene", "Sample", "Expression"))
}

RNAseq_heatmap_plotly_ignoregroup <- function(deseq2res, n = 100) {
  # Extract expression data
  # Assuming expression data starts from the 9th column to the second last column
  expr_data <- deseq2res[, 9:(ncol(deseq2res) - 1)]
  # Calculate variance for each gene across samples
  gene_variances <- apply(expr_data, 1, var, na.rm = TRUE)
  # Add variance as a new column to the results
  deseq2res$variance <- gene_variances
  # Select top 'n' genes with the highest variance
  top_genes <- deseq2res[order(deseq2res$variance, decreasing = TRUE), ][1:min(n, nrow(deseq2res)), ]
  # Define the heatmap title
  heatmap_title <- paste("Heatmap of Top", min(n, nrow(deseq2res)), "High Variance Genes")
  # Prepare the expression matrix for the heatmap
  # Adding 1 to avoid log10(0) issues and applying log10 transformation
  heatmap_matrix <- log10(top_genes[, 9:(ncol(top_genes) - 2)] + 1)
  # Remove row names to use custom labels
  rownames(heatmap_matrix) <- NULL
  # Generate the heatmap using heatmaply
  heatmaply::heatmaply(
    heatmap_matrix,
    main = heatmap_title,
    method = "plotly",
    labRow = top_genes$gene_name,  # Assuming 'gene_name' column exists
    xlab = "Samples",
    ylab = "Genes",
    showticklabels = c(TRUE, FALSE),
    show_dendrogram = c(FALSE, TRUE),
    key.title = "Scaled\nExpression\n(log10)",
    label_names = c("Gene", "Sample", "Expression"),
    fontsize_row = 10,               # Adjust as needed for readability
    fontsize_col = 10,               # Adjust as needed for readability
    colors = viridis::viridis(256)    # Optional: Use a color palette for better visualization
  )
}

RNAseq_PCA_ggplot2 <- function(deseq_output, samples2) {
  t <- t(deseq_output$norm_count[, c(-1,-2)])
  prin_comp <- prcomp(t, rank. = 2)
  components <- prin_comp[["x"]]
  components <- data.frame(components)
  components <- cbind(components, rownames(t))
  labels <- samples2$sample
  labels_group <- samples2$group
  if (length(unique(labels_group)) > 9) {
    library(Polychrome)
    mypalette <- as.vector(glasbey.colors(length(unique(labels_group))))
  } else {
    mypalette <- "Set1"
  }
  fig <- ggplot(components, aes(PC1, PC2, color=labels_group)) +
    geom_point(size=3) + theme_bw() + scale_fill_manual(values = mypalette) +
    labs(color = "Groups") + ggtitle("PCA") +
    theme(plot.title = element_text(hjust = 0.5))
  fig
}

RNAseq_volcano_ggplot2 <- function(deseq_output) {
  pal <- c("red", "royalblue")
  pal <- setNames(pal, c("Sig. genes", "Non-sig."))
  xmax <- max(abs(deseq_output$deseq2res$log2FoldChange)) * 1.1
  ymax <- max(-log10(deseq_output$deseq2res$padj)[is.finite(-log10(deseq_output$deseq2res$padj))]) * 1.1
  fig <- ggplot(deseq_output$deseq2res, aes(log2FoldChange, -log10(padj), color=sig)) +
        geom_point(size=1, alpha=0.2) + theme_bw() + scale_color_manual(values=pal) +
        labs(color = "Groups") + ggtitle("Volcano plot") +
        xlab("Fold Change (log2)") + ylab("adjusted p-value (-log10)") +
        theme(plot.title = element_text(hjust = 0.5)) +
        xlim(-xmax, xmax) + ylim(0, ymax)
  fig
}

RNAseq_maplot_ggplot2 <- function(deseq_output) {
  pal <- c("red", "royalblue")
  pal <- setNames(pal, c("Sig. genes", "Non-sig."))
  ymax <- max(abs(deseq_output$deseq2res$log2FoldChange)) * 1.01

  fig <- ggplot(deseq_output$deseq2res, aes(log2(baseMean), log2FoldChange, color=sig)) +
        geom_point(size=1, alpha=0.2) + theme_bw() + scale_color_manual(values=pal) +
        labs(color = "Groups") + ggtitle("MA plot") +
        xlab("Expresion Mean (log2)") + ylab("Fold Change (log2)") +
        theme(plot.title = element_text(hjust = 0.5)) +
        ylim(-ymax, ymax)
  fig
}

RNAseq_heatmap_ggplot2 <- function(deseq_output) {
  margin_spacer <- function(x) {
    # where x is the column in your dataset
    left_length <- nchar(levels(factor(x)))[1]
    if (left_length > 8) {
      return((left_length - 8) * 4)
    }
    else
      return(0)
  }
  sig_genes <- deseq_output$deseq2res[deseq_output$deseq2res$sig == "Sig. genes", ]
  heatmap_title <- "Heatmap of sig. DE genes"
  if (dim(sig_genes)[1] < 5) {
    sig_genes <- deseq_output$deseq2res[order(deseq_output$deseq2res$padj),][1:100,]
    heatmap_title <- "Heatmap of top 100 genes ranked by adj. p-value"
  }
  samples_names <- colnames(sig_genes)[9:(dim(sig_genes)[2]-1)]
  heatmap_t <- scale(log10(sig_genes[,9:(dim(sig_genes)[2]-1)]+1))
  ord <- hclust( dist(heatmap_t, method = "euclidean"), method = "ward.D" )$order

  heatmap_t <- cbind(sig_genes$gene_id, as.data.frame(heatmap_t))
  colnames(heatmap_t)[1] <- "gene_id"
  heatmap_t <- pivot_longer(heatmap_t, cols=2:(dim(heatmap_t)[2]), names_to="sample", values_to="Expression")
  heatmap_t$gene_id <- factor( heatmap_t$gene_id, levels = sig_genes$gene_id[ord])
  heatmap_t$sample <- factor( heatmap_t$sample, levels = samples_names)

  fig <- ggplot(heatmap_t, aes(sample, gene_id, fill=Expression)) +
        geom_tile() + ggtitle(heatmap_title) +
        ylab("Genes") + xlab("Samples") + scale_fill_viridis() +
        theme(plot.title = element_text(hjust = 0.5),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              plot.margin = margin(l = 0 + margin_spacer(heatmap_t$sample)))
  fig
}

RNAseq_heatmap_ggplot2_ignoregroup <- function(deseq_output, n = 100) {
  # Helper function to calculate left margin spacer based on sample name length
  margin_spacer <- function(x) {
    # x is the vector of sample names
    max_length <- max(nchar(as.character(x)), na.rm = TRUE)
    if (max_length > 8) {
      return((max_length - 8) * 4)
    } else {
      return(0)
    }
  }
  # Extract DESeq2 results
  deseq2res <- deseq_output$deseq2res
  # Check if 'gene_id' column exists
  if (!"gene_id" %in% colnames(deseq2res)) {
    stop("The 'deseq2res' data frame must contain a 'gene_id' column.")
  }
  # Extract expression data
  # Assuming expression data starts from the 9th column to the second last column
  expr_data <- deseq2res[, 9:(ncol(deseq2res) - 2)]
  # Calculate variance for each gene across samples
  gene_variances <- apply(expr_data, 1, var, na.rm = TRUE)
  # Add variance as a new column to the results
  deseq2res <- deseq2res %>%
    mutate(variance = gene_variances)
  # Select top 'n' genes with the highest variance
  top_genes <- deseq2res %>%
    arrange(desc(variance)) %>%
    slice_head(n = min(n, nrow(deseq2res)))
  # Define the heatmap title
  heatmap_title <- paste("Heatmap of Top", min(n, nrow(deseq2res)), "High Variance Genes")
  # Extract sample names
  samples_names <- colnames(top_genes)[9:(ncol(top_genes) - 1)]
  # Prepare the expression matrix
  # Adding 1 to avoid log10(0) issues and applying log10 transformation
  heatmap_t <- top_genes[, 9:(ncol(top_genes) - 1)] %>%
    log10(. + 1) %>%
    scale(center = TRUE, scale = TRUE)  # Standardize the data
  # Perform hierarchical clustering to order genes
  ord <- hclust(dist(heatmap_t, method = "euclidean"), method = "ward.D")$order
  # Combine gene_id with the scaled expression data
  heatmap_df <- top_genes %>%
    slice(ord) %>%
    select(gene_id) %>%
    bind_cols(as.data.frame(heatmap_t)) %>%
    pivot_longer(
      cols = -gene_id,
      names_to = "sample",
      values_to = "Expression"
    )
  # Set factor levels for proper ordering in the heatmap
  heatmap_df$gene_id <- factor(heatmap_df$gene_id, levels = top_genes$gene_id[ord])
  heatmap_df$sample <- factor(heatmap_df$sample, levels = samples_names)
  # Create the heatmap using ggplot2
  fig <- ggplot(heatmap_df, aes(x = sample, y = gene_id, fill = Expression)) +
    geom_tile() +
    ggtitle(heatmap_title) +
    ylab("Genes") +
    xlab("Samples") +
    scale_fill_viridis(name = "Scaled\nExpression\n(log10)", option = "viridis") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
      plot.margin = unit(c(1, 1, 1, margin_spacer(heatmap_df$sample)), "pt")
    )
  return(fig)
}

###########################################################
## Output tables
###########################################################

table_all_normalized_quantified_values <- function(normalized_counts) {
  datatable( normalized_counts ,
             extensions = c("Buttons" , "FixedColumns"),
             filter = 'top',
             options = list( autoWidth = TRUE ,
                             dom = 'Blftip',
                             pageLength = 10,
                             searchHighlight = FALSE,
                             buttons = c('copy', 'csv', 'print'),
                             scrollX = TRUE,
                             fixedColumns = list(leftColumns = 1)),
             class = c('compact cell-border stripe hover') ,
             rownames = FALSE) %>% formatRound(columns = c(-1,-2),digits=2)
}

table_diffexp_statistics <- function(deseq2res) {
  datatable( deseq2res ,
             extensions = c("Buttons" , "FixedColumns"),
             filter = 'top',
             options = list( autoWidth = TRUE ,
                             dom = 'Blftip',
                             pageLength = 10,
                             searchHighlight = FALSE,
                             buttons = c('copy', 'csv', 'print'),
                             scrollX = TRUE,
                             fixedColumns = list(leftColumns = 2)),
             class = c('compact cell-border stripe hover') ,
             rownames = FALSE) %>% formatRound(c(-1, -2), 2)
}

table_sig_genes <- function(res_sig, rowname=F) {
  res <- subset(res_sig, select = -c(sig) )
  datatable( res ,
             extensions = c("FixedColumns"),
             filter = 'top',
             options = list( autoWidth = TRUE ,
                             dom = 'Blftip',
                             pageLength = 10,
                             searchHighlight = FALSE,
                            #  buttons = c('copy', 'csv', 'print'),
                             scrollX = TRUE,
                             fixedColumns = list(leftColumns = 2)),
             class = c('compact cell-border stripe hover') ,
             rownames = rowname) %>% formatRound(c(-1, -2), 2)
}
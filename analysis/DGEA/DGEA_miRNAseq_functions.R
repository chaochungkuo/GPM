###########################################################
## miRNAseq
###########################################################

miRNAseq_deseq2 <- function(FILE_counts_mature, FILE_counts_hairpin, samples2, sizeFactorSpikeIn=FALSE, paired=FALSE) {
  df_counts1 <- t(read.csv(FILE_counts_mature, row.names = 1, header=T))[, samples2$sample]
  df_counts2 <- t(read.csv(FILE_counts_hairpin, row.names = 1, header=T))[, samples2$sample]
  # df_counts <- df_counts1
  df_counts <- rbind(df_counts1, df_counts2)
  if (paired) {
    dds <- DESeqDataSetFromMatrix(countData = df_counts, colData = samples2, design = ~batch + group)
  } else {
    dds <- DESeqDataSetFromMatrix(countData = df_counts, colData = samples2, design = ~group)
  }
  
  ## With QIAseq miRNA-Seq Spike-in
  if (sizeFactorSpikeIn) {
    sizeFactors(dds) <- sizeFactorSpikeIn
  } #################
  dds <- DESeq(dds)
  res <- as.data.frame(results(dds))
  res <- cbind(data.frame(gene_name=rownames(res)), res)
  res$sig <- "Non-sig."
  res$sig[res$padj < CUTOFF_ADJP] <- "Sig."
  ## With QIAseq miRNA-Seq Spike-in
  if (sizeFactorSpikeIn) {
    res$sig[grep("UniSP",res$gene_name)] <- "Spikein"
  } #################
  norm_counts <- counts(dds, normalized=TRUE)
  res_combined <- merge(res, norm_counts, by.x= "row.names", by.y="row.names", all.y=T)
  colnames(res_combined)[1] <- "gene_name"
  res_sig <- res[res$padj < CUTOFF_ADJP,]
  res_sig <- res_sig[complete.cases(res_sig), ]
  output <- list(res=res,
                 norm_counts=norm_counts,
                 res_combined=res_combined,
                 res_sig=res_sig)
  return(output)
}

miRNAseq_PCA_plotly <- function(normalized_counts, samples) {
  t <- t(normalized_counts)
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

miRNAseq_3D_PCA_plotly <- function(normalized_counts, samples) {
  t <- t(normalized_counts)
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

miRNAseq_volcano_plotly <- function(res_combined) {
  pal <- c("red", "royalblue")
  pal <- setNames(pal, c("Sig.", "Non-sig."))
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

miRNAseq_maplot_plotly <- function(res_combined) {
  pal <- c("red", "royalblue")
  pal <- setNames(pal, c("Sig.", "Non-sig."))
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

miRNAseq_heatmap_plotly <- function(res_combined) {
  res_reorder <- res_miRNA$res_combined[res_miRNA$res_combined$sig == "Sig.", ]
  heatmap_t <- log10(res_reorder[,10:(dim(res_reorder)[2])]+1)
  rownames(heatmap_t) <- c()
  heatmaply(heatmap_t, main = "Heatmap of DE miRNAs",
            method = "plotly", #labRow=res_reorder$gene_name,
            xlab = "Samples", ylab = "miRNA",
            showticklabels = c(TRUE, FALSE), show_dendrogram = c(FALSE, TRUE),
            key.title = "Scaled\nexpression\nin log10 scale",
            label_names = c("gene_name", "sample", "Expression"))
}

miRNAseq_PCA_ggplot2 <- function(res_miRNA, samples2){
  t <- t(res_miRNA$norm_counts)
  prin_comp <- prcomp(t, rank. = 2)
  components <- prin_comp[["x"]]
  components <- data.frame(components)
  components <- cbind(components, rownames(t))
  labels <- samples2$sample
  labels_group <- samples2$group

  fig <- ggplot(components, aes(PC1, PC2, color=labels_group)) +
        geom_point(size=3) + theme_bw() + scale_color_brewer(palette = "Set1") +
        labs(color = "Groups") + ggtitle("PCA") +
        theme(plot.title = element_text(hjust = 0.5))
  fig
}

miRNAseq_volcano_ggplot2 <- function(res_miRNA){
  pal <- c("red", "royalblue")
  pal <- setNames(pal, c("Sig.", "Non-sig."))
  xmax <- max(abs(res_miRNA$res$log2FoldChange)) * 1.1
  ymax <- max(-log10(res_miRNA$res$padj)[is.finite(-log10(res_miRNA$res$padj))]) * 1.1
  fig <- ggplot(res_miRNA$res, aes(log2FoldChange, -log10(padj), color=sig)) +
        geom_point(size=1, alpha=0.2) + theme_bw() + scale_color_manual(values=pal) +
        labs(color = "Groups") + ggtitle("Volcano plot") +
        xlab("Fold Change (log2)") + ylab("adjusted p-value (-log10)") +
        theme(plot.title = element_text(hjust = 0.5)) +
        xlim(-xmax, xmax) + ylim(0, ymax)
  fig
}

miRNAseq_maplot_ggplot2 <- function(res_miRNA){
  pal <- c("red", "royalblue")
  pal <- setNames(pal, c("Sig.", "Non-sig."))
  ymax <- max(abs(res_miRNA$res$log2FoldChange)) * 1.01

  fig <- ggplot(res_miRNA$res, aes(log2(baseMean), log2FoldChange, color=sig)) +
        geom_point(size=1, alpha=0.2) + theme_bw() + scale_color_manual(values=pal) +
        labs(color = "Groups") + ggtitle("MA plot") +
        xlab("Expresion Mean (log2)") + ylab("Fold Change (log2)") +
        theme(plot.title = element_text(hjust = 0.5)) +
        ylim(-ymax, ymax)
  fig
}

miRNAseq_heatmap_ggplot2 <- function(res_miRNA){
  margin_spacer <- function(x) {
    # where x is the column in your dataset
    left_length <- nchar(levels(factor(x)))[1]
    if (left_length > 8) {
      return((left_length - 8) * 4)
    }
    else
      return(0)
  }
  res_reorder <- res_miRNA$res_combined[res_miRNA$res_combined$sig == "Sig.", ]
  samples_names <- colnames(res_reorder)[10:(dim(res_reorder)[2])]
  heatmap_t <- scale(log10(res_reorder[,10:(dim(res_reorder)[2])]+1))
  ord <- hclust( dist(heatmap_t, method = "euclidean"), method = "ward.D" )$order

  heatmap_t <- cbind(res_reorder$gene_name, as.data.frame(heatmap_t))
  colnames(heatmap_t)[1] <- "gene_name"
  heatmap_t <- pivot_longer(heatmap_t, cols=2:(dim(heatmap_t)[2]), names_to="sample", values_to="Expression")
  heatmap_t$gene_name <- factor( heatmap_t$gene_name, levels = res_reorder$gene_name[ord])
  heatmap_t$sample <- factor( heatmap_t$sample, levels = samples_names)

  fig <- ggplot(heatmap_t, aes(sample, gene_name, fill=Expression)) +
        geom_tile() + ggtitle("Heatmap of DE genes") +
        ylab("Genes") + xlab("Samples") + scale_fill_viridis() +
        theme(plot.title = element_text(hjust = 0.5),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              plot.margin = margin(l = 0 + margin_spacer(heatmap_t$sample)))
  fig
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
             rownames = FALSE) %>% formatRound(columns = c(-1),digits=2)
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
                             fixedColumns = list(leftColumns = 1)),
             class = c('compact cell-border stripe hover') ,
             rownames = FALSE) %>% formatRound(c(-1), 2)
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
                             fixedColumns = list(leftColumns = 1)),
             class = c('compact cell-border stripe hover') ,
             rownames = rowname) %>% formatRound(c(-1), 2)
}
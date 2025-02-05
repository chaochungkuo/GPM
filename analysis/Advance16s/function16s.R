###########################################################
## Collection of functions during 16s analysis
###########################################################
library(plotly)
library(heatmaply)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(RColorBrewer)
library(lefser)
library(qiime2R)
library(phyloseq)
library(microbiome)
library(microbiomeMarker)
library(purrr)
library(dplyr)
library(forcats)
library(microViz)

##function to load result from nf-core include species level
load_physeq_full <- function(ASVs, meta, taxonomy, tree)
{
  
  tax_table <- do.call(rbind, strsplit(as.character(taxonomy$data$Taxon), ";"))
  tax_table <- tax_table[,1:7]
  
  colnames(tax_table) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  rownames(tax_table) <- taxonomy$data$Feature.ID
  
  
  ##change species terms according to scientific name guide
  tax_table[,7] <- apply( tax_table[,6:7] , 1 , paste , collapse = " " )
  
  physeq <- phyloseq::phyloseq(
    otu_table(ASVs$data, taxa_are_rows = T),
    phy_tree(tree$data),
    tax_table(tax_table),
    sample_data(meta)
  )
  return(physeq)
  
}

## clean relative abundant table for positive control sample for plot 
#species with confidence over 0.95 in the positive control sample are presented 

pos_bartab(){
    bar.tab <- read.delim(paste0(DIR_result,"/qiime2/rel_abundance_tables/rel-table-ASV_with-DADA2-tax.tsv"))
    bar.tab <- bar.tab[,c(-1,-9)] # remove some annotations
    bar.tab <- bar.tab[(bar.tab[,8] >= 0.95),]
    bar.tab$species <- apply(bar.tab[,6:7],1, paste,  collapse=".")
    bar.tab <- bar.tab[,c(122,11)] #point the column of positive control
    bar.tab <- aggregate(bar.tab[,-1],bar.tab["species"],sum)
    return(bar.tab)
}

## function to plot relative abundancy of positive control
rel_abundance_bar_plotly <- function(rel_abund) {
  
  rel_abund <- rel_abund[!(rel_abund[,-1] == 0),]
  #rel_abund[,1] <- str_extract(rel_abund[,1],"\\w+$")
  rel_abund <- setNames(data.frame(t(rel_abund[ , -1])), rel_abund[ ,1])
  
  fig <- data.table::melt(as.matrix(rel_abund)) %>%
    plot_ly(x = ~Var1, y = ~value, type = 'bar', 
            name = ~Var2, color = ~Var2, colors= viridis(n = 256, alpha = 1, begin = 1,end = 0, option = "viridis")) %>%
    layout(yaxis = list(title = 'Rel_abundance'), xaxis = list( title = 'Positive control'), barmode = 'stack')
  
  fig
}

 
## function to run linear discriminal analysis
run_cus_lefse <- function(physeq, groupN, subgroupN){
  #filter the taxa 
  physeq.filter <- filter_taxa(physeq, function(x)  sum(x> 10) >0, TRUE)
  
  mm_lefse<- run_lefse(
    physeq.filter,
    group = groupN,
    subgroup = subgroupN,
    norm = "CPM",
    multigrp_strat = T,
    only_same_subgrp = T,
    lda_cutoff = 0,
    kw_cutoff = 1,
    wilcoxon_cutoff = 1
    #multigrp_strat = T
  )
  
  return(mm_lefse)
}

## visualization functions ##
#############################


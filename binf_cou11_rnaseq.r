# File      binf_cou11_rnaseq
# Version   0.9
# Date      23/05/2018
# Authors   Sjors Bongers, Daan Gilissen, Martijn Landman, Koen Rademaker, Ronald van den Hurk

source('http://bioconductor.org/biocLite.R')
# (Optional) installation of packages.
# biocLite('edgeR')
# biocLite('xlsx')
# biocLite('KEGGREST')
library('edgeR')
library('xlsx')
library(KEGGREST)

Counts <- read.delim('Data/RNA-Seq-counts.txt', header=TRUE, skip=1, row.names=1)
Annotation <- read.delim('Data/RNA-Seq-annotation.txt', header=TRUE, skip=1, row.names=1)
CPM = 10
PCH_1 = 21
PCH_2 = 23

CreateGroup <- function(conditions){
  # Store experimental conditions.
  exp <- rep(conditions, each=2)
  group <- factor(exp)
  
  return(group)
}

CreateModel <- function(strain, data, group){
  # Group samples by condition into design matrix.
  design <- model.matrix(~0+group, data=data$samples)
  colnames(design) <- levels(data$samples$group)
  
  # Create model for top differentially expressed genes.
  fit <- glmFit(data, design)
  if(strain == 'WCFS1'){
    mc <- makeContrasts(exp.r=WCFS1.glc-WCFS1.rib, levels=design)
  } else if(strain == 'NC8'){
    mc <- makeContrasts(exp.r=NC8.glc-NC8.rib, levels=design)
  }
  fit <- glmLRT(fit, contrast=mc)
  
  return(fit)
}

DataProcessing <- function(group, start, stop, cpm_filter){
  # Create DGEList object for storage of RNA-Seq data.
  y <- DGEList(counts=Counts[,start:stop], group=group)
  
  # Filter out genes below cpm_filter.
  keep.genes <- rowSums(cpm(y)>cpm_filter) >= 2
  y <- y[keep.genes,]
  y$samples$lib.size <- colSums(y$counts)
  
  # Determine scale factors using Trimmed Mean of M-values (TMM).
  y <- calcNormFactors(y, method='TMM')
  
  # Group samples by condition into design matrix.
  design <- model.matrix(~0+group, data=y$samples)
  colnames(design) <- levels(y$samples$group)
  
  # Estimate dispersions.
  y <- estimateGLMCommonDisp(y,design)
  y <- estimateGLMTrendedDisp(y,design, method='power')
  y <- estimateGLMTagwiseDisp(y,design)
  
  return(y)
}

PlotSampleDistances <- function(title, data, group){
  # Set up colors and symbols for plot.
  if(length(levels(group)) == 4){
    colors <- rep(c('red', 'red', 'blue', 'blue'), 2)
  } else if(length(levels(group)) == 2){
    colors <- rep(c('red', 'blue'), 2)
  }
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  pch <- rep(c(PCH_1, PCH_2), length(levels(group)))
  pch_legend <- rep(c(16, 18), length(levels(group))/2)
  
  # Visualize plot.
  plotMDS(data, bg=colors[group], cex=2, col=1, pch=pch[group], xlab='Dimension 1', ylab='Dimension 2')
  title(title, line=0.5)
  legend('topright', col=colors, inset=c(-0.33,0), legend=levels(group), ncol=1, pch=pch_legend, title='Samples')
}

GetPathwaysForGenes <- function(genes){
  # Set up dataframe.
  rows_genes <- rownames(genes)
  n_pathways <- 0
  for(i in 1:length(rows_genes)){
    gene <- rows_genes[i]
    try(query <- keggGet(c(paste('lpl:', gene, sep=''))), silent=F)
    if(exists('query')){
      pathways <- query[[1]]$PATHWAY
      if(!is.null(pathways)){
        if(length(pathways) > n_pathways){
          n_pathways = length(pathways)
        }
      }
    }
  }
  pathways_genes <- data.frame(matrix(ncol = length(rows_genes), nrow = n_pathways))
  colnames(pathways_genes) <- rows_genes
  
  # Store pathways per gene in dataframe.
  for(i in 1:length(rows_genes)){
    gene <- rows_genes[i]
    try(query <- keggGet(c(paste('lpl:', gene, sep=''))), silent=F)
    if(exists('query')){
      pathways <- query[[1]]$PATHWAY
      if(!is.null(pathways)){
        for(j in 1:length(pathways)){
          pathways_genes[j, i] = pathways[j]
        }
      }
    }
  }
  
  return(pathways_genes)
}

DetermineDEGenes <- function(fit, n_results){
  # Determine top differentially expressed genes.
  top_DE_genes <- topTags(fit, n=n_results)
  return(top_DE_genes)
}

DeterminePathwayOverrep <- function(fit, n_results){
  # Determine overrepresentation of genes in KEGG pathways.
  kegg_pathways <- kegga(fit, species.KEGG='lpl')
  top_OR_pathways <- topKEGG(kegg_pathways, number=n_results)
  
  return(top_OR_pathways)
}

AnnotateDEGEnes <- function(genes){
  # Annotate DE genes.
  annotated_data <- cbind(genes, Annotation[rownames(genes),])
  return(annotated_data)
}

WriteResults <- function(file_name, annotated_results, sheet_name_1, or_pathways, sheet_name_2, pathways_de_genes, sheet_name_3){
  # Write results to sheets in Exel file.
  write.xlsx(annotated_results, file=file_name, sheetName=sheet_name_1, col.names=TRUE, row.names=TRUE, append=TRUE, showNA=TRUE)
  write.xlsx(or_pathways, file=file_name, sheetName=sheet_name_2, col.names=TRUE, row.names=TRUE, append=TRUE, showNA=TRUE)
  write.xlsx(pathways_de_genes, file=file_name, sheetName=sheet_name_3, col.names=TRUE, row.names=FALSE, append=TRUE, showNA=FALSE)
}

# Visualize and check separation of samples on growth medium and strain.
All_group <- CreateGroup(c('WCFS1.glc', 'WCFS1.rib', 'NC8.glc', 'NC8.rib'))
All_data <- DataProcessing(All_group, 1, 8, CPM)
PlotSampleDistances('Distances between RNA-Seq samples', All_data, All_group)

# Run analysis for WCFS1.
WCFS1_group <- CreateGroup(c('WCFS1.glc', 'WCFS1.rib'))
WCFS1_data <- DataProcessing(WCFS1_group, 1, 4, CPM)
fit <- CreateModel('WCFS1', WCFS1_data, WCFS1_group)
de_genes <- DetermineDEGenes(fit, nrow(WCFS1_data))
pathways_de_genes <- GetPathwaysForGenes(de_genes)
overrep_pathways <- DeterminePathwayOverrep(fit, Inf)
annotated_results <- AnnotateDEGEnes(de_genes)
PlotSampleDistances('Distances between WCFS1 RNA-Seq samples', WCFS1_data, WCFS1_group)
WriteResults('Results/RNA_Seq_analysis_results.xlsx', annotated_results, 'WCFS1 DE genes', overrep_pathways, 'WCFS1 Overrep pathways', pathways_de_genes, 'WCFS1 DE genes pathways')

# Run analysis for NC8.
NC8_group <- CreateGroup(c('NC8.glc', 'NC8.rib'))
NC8_data <- DataProcessing(NC8_group, 5, 8, CPM)
fit <- CreateModel('NC8', NC8_data, NC8_group)
de_genes <- DetermineDEGenes(fit, nrow(NC8_data))
pathways_de_genes <- GetPathwaysForGenes(de_genes)
overrep_pathways <- DeterminePathwayOverrep(fit, Inf)
annotated_results <- AnnotateDEGEnes(de_genes)
PlotSampleDistances('Distances between NC8 RNA-Seq samples', NC8_data, NC8_data)
WriteResults('Results/RNA_Seq_analysis_results.xlsx', annotated_results, 'NC8 DE genes', overrep_pathways, 'NC8 Overrep pathways', pathways_de_genes, 'NC8 DE genes pathways')

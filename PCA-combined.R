# File      binf_cou11_rnaseq
# Version   0.8
# Date      17/05/2018
# Authors   Sjors Bongers, Daan Gilissen, Martijn Landman, Koen Rademaker, Ronald van den Hurk

source('http://bioconductor.org/biocLite.R')
# (Optional) installation of packages.
# biocLite('edgeR')
# biocLite('xlsx')
# biocLite('KEGGREST')
library('edgeR')
library('xlsx')
library(KEGGREST)

#Links aangepast aan mijn PC (want windows)
Counts <- read.delim('C:\\Users\\Sjors\\Documents\\Github\\RNA-seq-analysis\\Data\\RNA-Seq-counts.txt', header=TRUE, skip=1, row.names=1)
Annotation <- read.delim('C:\\Users\\Sjors\\Documents\\Github\\RNA-seq-analysis\\Data\\RNA-Seq-annotation.txt', header=TRUE, skip=1, row.names=1)
CPM = 10


exp <- c("WCFS1.glc","WCFS1.glc","WCFS1.rib","WCFS1.rib", "NC8.glc", "NC8.glc", "NC8.rib", "NC8.rib")
group <- factor(exp)

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

norm_data <- DataProcessing(group, 1, 8, CPM)
plotMDS(y)
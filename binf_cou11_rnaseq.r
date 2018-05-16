# File      binf_cou11_rnaseq
# Version   0.7
# Date      14/05/2018
# Authors   Sjors Bongers, Daan Gilissen, Martijn Landman, Koen Rademaker, Ronald van den Hurk

source('http://bioconductor.org/biocLite.R')
# (Optional) installation of packages.
# biocLite('edgeR')
# biocLite('xlsx')
# biocLite('KEGGREST')
library('edgeR')
library('xlsx')
library(KEGGREST)

counts <- read.delim('Data/RNA-Seq-counts.txt', header=TRUE, skip=1, row.names=1)
annotation <- read.delim('Data/RNA-Seq-annotation.txt', header=TRUE, skip=1, row.names=1)
cpm = 10

RNA_seq_analysis <- function(exp_val_1, exp_val_2, index_1, index_2, cpm_filter, pdf_name, save_pdf, xlsx_name, data_sheet_name, kegg_sheet_name, pathways_sheet_name, save_sheet){
  # Create DGEList object for storing data.
  exp <-c(exp_val_1, exp_val_1, exp_val_2, exp_val_2)
  group <- factor(exp)
  y <- DGEList(counts=counts[,index_1:index_2], group=group)
  
  # Select genes with a CPM above or equal to cpm_filter.
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
  
  # Plot normalized data.
  if(save_pdf){
    pdf(pdf_name)
    pch <- c(15, 17, 15, 17)
    colors <- rep(c('red', 'blue'), 2)
    plotMDS(y, col=colors[group], pch=pch[group])
    legend('topright', legend=levels(group), pch=pch, col=colors, ncol=1)
    plotBCV(y)
    dev.off()
  }
  
  # Determine differentially expressed genes.
  fit <- glmFit(y,design)
  if(exp_val_1 == 'WCFS1.glc' & exp_val_2 == 'WCFS1.rib'){
    mc <- makeContrasts(exp.r=WCFS1.glc-WCFS1.rib, levels=design)
  }
  if(exp_val_1 == 'NC8.glc' & exp_val_2 == 'NC8.rib'){
    mc <- makeContrasts(exp.r=NC8.glc-NC8.rib, levels=design)
  }
  fit <- glmLRT(fit, contrast=mc)
  res_toptags <- topTags(fit, n=nrow(y))
  
  # Set up dataframe for top 10 differentially expressed genes and corresponding pathways.
  res_top10tags <- topTags(fit)
  rows_top10tags <- rownames(res_top10tags)
  max_n_pathways = 0
  for(i in 1:length(rows_top10tags)){
    gene <- rows_top10tags[i]
    query <- keggGet(c(paste('lpl:', gene, sep='')))
    pathways <- query[[1]]$PATHWAY
    if(!is.null(pathways)){
      if(length(pathways) > max_n_pathways){
        max_n_pathways = length(pathways)
      }
    }
  }
  res_pathways <- data.frame(matrix(ncol = length(rows_top10tags), nrow = max_n_pathways))
  colnames(res_pathways) <- rows_top10tags
  
  # Store pathways for top 10 differentially expressed genes.
  for(i in 1:length(rows_top10tags)){
    gene <- rows_top10tags[i]
    query <- keggGet(c(paste('lpl:', gene, sep='')))
    pathways <- query[[1]]$PATHWAY
    if(!is.null(pathways)){
      for(j in 1:length(pathways)){
        res_pathways[j, i] = pathways[j]
      }
    }
  }
  
  # Combine top differentially expressed genes with annotation.
  combined_data <- cbind(res_toptags, annotation[rownames(res_toptags),])
  
  # KEGG pathway analysis.
  kegg <- kegga(fit, species.KEGG='lpl')
  res_kegg <- topKEGG(kegg)
  
  # Save data to Excel file.
  if(save_sheet){
    write.xlsx(combined_data, file=xlsx_name, sheetName=data_sheet_name, col.names=TRUE, row.names=TRUE, append=TRUE, showNA=TRUE)
    write.xlsx(res_kegg, file=xlsx_name, sheetName=kegg_sheet_name, col.names=TRUE, row.names=TRUE, append=TRUE, showNA=TRUE)
    write.xlsx(res_pathways, file=xlsx_name, sheetName=pathways_sheet_name, col.names=TRUE, row.names=FALSE, append=TRUE, showNA=FALSE)
  }
}

WCFS1 <- RNA_seq_analysis('WCFS1.glc', 'WCFS1.rib', 1, 4, cpm, 'Results/WCFS1 results.pdf', FALSE, 'Results/Annotated DE genes and pathways.xlsx', 'WCFS1 Top D.E. genes', 'WCFS1 Top pathways', 'WCFS1 Top 10 genes pathways', FALSE)
NC8 <- RNA_seq_analysis('NC8.glc', 'NC8.rib', 5, 8, cpm, 'Results/NC8 results.pdf', FALSE, 'Results/Annotated DE genes and pathways.xlsx', 'NC8 Top D.E. genes', 'NC8 Top pathways', 'NC8 Top 10 genes pathways', FALSE)

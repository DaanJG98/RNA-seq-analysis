# Version   0.14
# Date      04/06/2018

# INSTALLATION (OPTIONAL) AND LOADING REQUIRED PACKAGES
source('http://bioconductor.org/biocLite.R')
if(is.element('biocLite', installed.packages()[,1])==FALSE){
  biocLite()
}
required_libraries <- c('DBI', 'edgeR', 'gplots', 'KEGGREST', 'KEGG.db', 'RColorBrewer', 'xlsx')
for(library in required_libraries){
  if(is.element(library, installed.packages()[,1])==FALSE){
    biocLite(pkgs=library)
  }
}
library('DBI')
library('edgeR')
library('gplots')
library('KEGGREST')
library('KEGG.db')
library('RColorBrewer')
library('xlsx')

# LOAD DATA
Counts <- read.delim('Data/RNA-Seq-counts.txt', header=TRUE, skip=1, row.names=1)
Annotation <- read.delim('Data/RNA-Seq-annotation.txt', header=TRUE, skip=1, row.names=1)

# SET GLOBAL VARIABLES
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

PlotHeatMap <- function(WCFS1_df, NC8_df, n_genes){
  # Cuts data frames down to selected number of genes.
  WCFS1_df <- WCFS1_df[1:n_genes,]
  NC8_df <- NC8_df[1:n_genes,]
  shared_genes <- intersect(row.names(WCFS1_df), row.names(NC8_df))
  # Reduces data frames to only contain shared genes.
  WCFS1_df <- as.data.frame(WCFS1_df[row.names(WCFS1_df) %in% shared_genes,])
  NC8_df <- as.data.frame(NC8_df[row.names(NC8_df) %in% shared_genes,])
  # Combines data frames into single data frame and provides annotation.
  annotation_df <- as.data.frame(Annotation[row.names(Annotation) %in% shared_genes,])
  annotation_df <- subset(annotation_df, select=c('name', 'subclass'))
  df <- cbind(WCFS1_df['logFC'], NC8_df['logFC'], row.names=shared_genes)
  annotation_df <- cbind(df, annotation_df[rownames(df),], row.names=shared_genes)
  # Visualizes heat map and sets column and row names.
  rownames(df) <- paste(annotation_df[,"name"], annotation_df[,"subclass"], sep=' - ')
  colnames(df) <- c('WCFS1', 'NC8')
  color_palette <- colorRampPalette(c("green","blue"))(n = 64)
  heatmap.2(as.matrix(df),
            adjCol=c(NA, 0.5),
            cexCol=1.5,
            col=color_palette,
            Colv=F,
            dendrogram='none',
            density.info='none',
            margins=c(3,22),
            notecol='black',
            srtCol=0,
            trace='none')
}

GetPathwaysForGenes <- function(genes){
  # Set up data frame.
  cols <- rownames(genes)
  n_pathways <- dbGetQuery(KEGG_dbconn(), 'SELECT COUNT(*) FROM pathway2name')[1,1]
  df <- data.frame(matrix(ncol=length(cols), nrow=n_pathways))
  colnames(df) <- cols
  # Store pathways per gene in data frame.
  max_n_pathways = 0
  for(index in 1:length(cols)){
    gene <- cols[index]
    try(query <- keggGet(c(paste('lpl:', gene, sep=''))), silent=F)
    if(exists('query')){
      pathways <- query[[1]]$PATHWAY
      if(!is.null(pathways)){
        for(index_2 in 1:length(pathways)){
          df[index_2, index] = pathways[index_2]
        }
        if(length(pathways) > max_n_pathways){
          max_n_pathways = length(pathways)
        }
      }
    }
  }
  df <- df[-c(max_n_pathways+1:nrow(df)), ]
  return(df)
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

AnnotateDEGenes <- function(genes){
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


# VISUALIZE AND CHECK SEPARATION OF SAMPLES ON GROWTH MEDIUM AND STRAIN
All_group <- CreateGroup(c('WCFS1.glc', 'WCFS1.rib', 'NC8.glc', 'NC8.rib'))
All_data <- DataProcessing(All_group, 1, 8, CPM)
PlotSampleDistances('Distances between RNA-Seq samples', All_data, All_group)

# PROCESS DATA - WCFS1
WCFS1_group <- CreateGroup(c('WCFS1.glc', 'WCFS1.rib'))
WCFS1_data <- DataProcessing(WCFS1_group, 1, 4, CPM)
WCFS1_fit <- CreateModel('WCFS1', WCFS1_data, WCFS1_group)
# PROCESS DATA - NC8
NC8_group <- CreateGroup(c('NC8.glc', 'NC8.rib'))
NC8_data <- DataProcessing(NC8_group, 5, 8, CPM)
NC8_fit <- CreateModel('NC8', NC8_data, NC8_group)

# DETERMINE DE GENES AND PATHWAYS - WCFS1
WCFS1_de_genes <- DetermineDEGenes(WCFS1_fit, nrow(WCFS1_data))
WCFS1_pathways_de_genes <- GetPathwaysForGenes(WCFS1_de_genes)
# DETERMINE DE GENES AND PATHWAYS - NC8
NC8_de_genes <- DetermineDEGenes(NC8_fit, nrow(NC8_data))
NC8_pathways_de_genes <- GetPathwaysForGenes(NC8_de_genes)

# VALIDATE RESULTS - WCFS1
WCFS1_overrep_pathways <- DeterminePathwayOverrep(WCFS1_fit, Inf)
WCFS1_annotated_results <- AnnotateDEGenes(WCFS1_de_genes)
# VALIDATE RESULTS - NC8
NC8_overrep_pathways <- DeterminePathwayOverrep(NC8_fit, Inf)
NC8_annotated_results <- AnnotateDEGenes(NC8_de_genes)

# PLOT DATA
PlotSampleDistances('Distances between WCFS1 RNA-Seq samples', WCFS1_data, WCFS1_group)
PlotSampleDistances('Distances between NC8 RNA-Seq samples', NC8_data, NC8_data)
PlotHeatMap(as.data.frame(WCFS1_de_genes$table), as.data.frame(NC8_de_genes$table), 50)

# EXPORT RESULTS
WriteResults('Results/RNA_Seq_analysis_results.xlsx', WCFS1_annotated_results, 'WCFS1 DE genes', WCFS1_overrep_pathways, 'WCFS1 Overrep pathways', WCFS1_pathways_de_genes, 'WCFS1 DE genes pathways')
WriteResults('Results/RNA_Seq_analysis_results.xlsx', NC8_annotated_results, 'NC8 DE genes', NC8_overrep_pathways, 'NC8 Overrep pathways', pathways_de_genes, 'NC8 DE genes pathways')

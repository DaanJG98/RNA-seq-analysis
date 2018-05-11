# File      binf_cou11_kegg
# Version   0.2
# Date      10/05/2018
# Authors   Sjors Bongers, Daan Gilissen, Martijn Landman, Koen Rademaker, Ronald van den Hurk

source('http://bioconductor.org/biocLite.R')
# (Optional) Installation of package.
#biocLite('KEGGREST')
library(KEGGREST)

# Create query
gene <- 'lp_0001'
query <- keggGet(c(paste('lpl:', gene, sep='')))

# Get gene name
query[[1]]$NAME
# Get pathway(s)
query[[1]]$PATHWAY
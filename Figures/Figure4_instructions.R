# Script to reproduce Figure 4 #
################################

# Load R 4.2.1
# Install singleCellTK 2.8.1

library(singleCellTK)
sce <- readRDS("tutorial_pbmc3k_qc.rds")
assay(sce, "counts") <- assay(sce, "decontXcounts")
reportSeurat(inSCE = sce, 
             biological.group = "decontX_clusters", 
             selected.markers = c("CD3E", "CD79A"))

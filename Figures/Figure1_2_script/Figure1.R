# Script to generate all the subplot for Figure 1 ####

# setup ####
library(singleCellTK)
library(celda)
library(ggplot2)
library(cowplot)

source("helper.R")
packageVersion("singleCellTK") # '2.8.1' at Jan 2nd, 2023 (scanpy branch)
seed <- 12345

# Import data ####
sce <- importExampleData("pbmc3k")


# Run QC ####
sce <- runCellQC(sce, algorithms = c("QCMetrics", "scDblFinder"),
                 sample = sce$Sample, mitoRef = "human", mitoIDType = "symbol",
                 mitoGeneLocation = "rownames", seed = seed)
sce <- runQuickUMAP(sce, initialDims = 10)

# Plot QC ####
## Plot per cell Generic QC ####
plots.percellQC <- plotRunPerCellQCResults(sce, combinePlot = "none",
                                           baseSize = 16)

pdf("percellQC.pdf", width = 10, height = 9)
plots.percellQC$Sum +
  theme(text = element_text(size = 28))
dev.off()
## Plot scDblFinder ####
plots.scDblFinder <- plotScDblFinderResults(sce, baseSize = 14, combinePlot = "none")
pdf("scDblFinder.pdf", width = 10, height = 9)
plots.scDblFinder$scatter_doubletCall
dev.off()



# A la carte Workflow ####
sce <- subsetSCECols(sce, colData = c("total > 600", "detected > 300",
                                      "mito_percent < 5"))
sce <- runNormalization(sce, useAssay = "counts", outAssayName = "logcounts",
                        normalizationMethod = "logNormCounts")
sce <- runModelGeneVar(sce, useAssay = "logcounts")
sce <- setTopHVG(sce, method = "modelGeneVar", hvgNumber = 2000,
                 featureSubsetName = "hvg2000")
sce <- scaterPCA(sce, useFeatureSubset = "hvg2000", seed = seed)
sce <- runUMAP(sce, useReducedDim = "PCA", initialDims = 10, seed = seed)
sce <- runScranSNN(sce, "PCA", clusterName = "cluster", nComp = 10,
                   weightType = "jaccard", k = 14, seed = seed)

# Plot A La Carte ####
## plot HVG ####
pdf("hvg.pdf", width = 10, height = 9)
plotTopHVG(sce, method = "modelGeneVar", hvgNumber = 2000, labelsCount = 5,
           textSize = 18, labelSize = 6)
dev.off()
## plot PCA Elbow ####
seuratObj <- convertSCEToSeurat(sce)
sce@metadata$seurat$obj <- seuratObj
reDim <- reducedDim(sce, "PCA")
colnames(reDim) <- paste0("PCA_", seq_len(length(colnames(reDim))))
rownames(reDim) <- gsub('_', '-', rownames(reDim))
sce@metadata$seurat$obj@reductions$pca <-
    Seurat::CreateDimReducObject(embeddings = reDim,
                                 key = "PCA_", assay = "RNA",
                                 loadings = attr(reDim, "rotation"),
                                 stdev = as.numeric(attr(reDim, "percentVar")))
pdf("elbowPlot.pdf", width = 10, height = 9)
plotSeuratElbow(sce, interactive = FALSE) +
    theme(text = element_text(size = 16),
          axis.text = element_text(size = 16))
dev.off()
## plot UMAP ####
pdf("UMAP.pdf", width = 10, height = 9)
plotSCEDimReduceColData(sce, "cluster", "UMAP", labelClusters = FALSE,
                        baseSize = 12) +
    theme(legend.position = "none")
dev.off()



# Seurat Curated Workflow ####
sce <- runSeuratNormalizeData(sce, useAssay = "counts")
sce <- runSeuratFindHVG(sce, method = "vst", hvgNumber = 2000,
                        useAssay = "counts")
sce <- setTopHVG(sce, method = "vst", hvgNumber = 2000,
                 featureSubsetName = "vst2000")
sce <- runSeuratPCA(sce, useAssay = "seuratNormData",
                    useFeatureSubset = "vst2000", scale = TRUE, seed = seed)
sce <- runSeuratFindClusters(sce, useAssay = "seuratNormData",
                             useReduction = "pca")
# Cluster variable name: "Seurat_louvain_Resolution0.8"
sce <- runSeuratFindMarkers(sce, allGroup = "Seurat_louvain_Resolution0.8")

# Plot Seurat Curated Workflow ####
## PCA scatter####
pdf("Seurat_PCA.pdf", width = 10, height = 9)
plotSeuratReduction(sce, useReduction = "pca") +
    theme(text = element_text(size = 28),
          axis.text = element_text(size = 28))
dev.off()
## PCA Loading Heatmap####
sce@metadata$seurat$heatmap_dimRed <- computeHeatmap(
    sce, useAssay = "seuratNormData", dims = 1:4, nfeatures = 16,
    reduction = "pca"
)
pdf("Seurat_PCA_Heatmap.pdf", width = 10, height = 8)
plot_grid(plotlist = sce@metadata$seurat$heatmap_dimRed)
dev.off()
## Seurat marker Ridge Plot ####
ridge_list <- plotSeuratGenes(
    sce, plotType = "ridge",
    features = sce@metadata$seuratMarkers$gene.id[1:4],
    groupVariable = "Seurat_louvain_Resolution0.8",
    ncol = 2, combine = FALSE
)
for (i in 1:4) {
    ridge_list[[i]] <- ridge_list[[i]] +
        theme(axis.text = element_text(size = 22),
              axis.title = element_text(size = 24),
              plot.title = element_text(size = 28),
              legend.position = "none")
}
svg("Seurat_Marker_RidgePlot.svg", width = 10, height = 8)
plot_grid(plotlist = ridge_list)
dev.off()



# Celda Curated Workflow Analysis ####
useAssay <- "counts"
altExpName <- "featureSubset"
sce <- setTopHVG(sce, method = "modelGeneVar", hvgNumber = 5000,
                 featureSubsetName = "hvg5000")
altExp(sce, altExpName) <- sce[rowSubset(sce, "hvg5000"), ]
moduleSplit <- recursiveSplitModule(sce, useAssay = useAssay,
                                    altExpName = altExpName, initialL = 10,
                                    maxL = 150, seed = seed)
L <- 80
temp <- subsetCeldaList(moduleSplit, list(L = L))
sce <- recursiveSplitCell(sce, useAssay = useAssay, altExpName = altExpName,
                          initialK = 3, maxK = 25, yInit = celdaModules(temp),
                          seed = seed)
K <- 14
sce <- subsetCeldaList(sce, list(L = L, K = K))
mod <- featureModuleLookup(sce, feature = c("CD3E"))
sce <- celdaUmap(sce, useAssay = useAssay, altExpName = altExpName)

# Plot Celda Analysis ####
## Module RPC Plot ####
pdf("celda_Module_RPC.pdf", width = 10, height = 6)
plotRPC(moduleSplit, altExpName = altExpName, sep = 10) +
    theme(text = element_text(size = 28)) +
    ylab("Rate of Perplexity Change")
dev.off()
## Cell Population RPC Plot ####
pdf("celda_cell_RPC.pdf", width = 10, height = 6)
plotRPC(sce, sep = 3) +
    theme(text = element_text(size = 28)) +
  ylab("Rate of Perplexity Change")
dev.off()
## Module Heatmap ####
pdf("celda_module_heatmap.pdf", width = 10, height = 6)
moduleHeatmap(sce, featureModule = mod, useAssay = useAssay,
              altExpName = altExpName, rowFontSize = 16, moduleLabelSize = 28)
dev.off()
## Module Expression UMAP ####
pdf("celda_module_UMAP.pdf", width = 10, height = 6)
plotDimReduceModule(sce, modules = mod, useAssay = useAssay,
                    altExpName = altExpName, reducedDimName = "celda_UMAP") +
    theme(text = element_text(size = 28))
dev.off()


# Scanpy Curated Workflow ####
set.seed(12345)
sce <- runScanpyNormalizeData(sce)
sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData")
sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData", use_highly_variable = TRUE)
sce <- runScanpyUMAP(sce, useReducedDim = "scanpyPCA", dims = 10)
sce <- runScanpyFindClusters(sce, dims = 10L, nNeighbors = 30, resolution = 0.8, method = "louvain")
sce <- runScanpyFindMarkers(sce, colDataName = "Scanpy_louvain_0.8")

pdf("scanpy_pca_loading.pdf", width = 10, height = 6)
plotScanpyPCLoading(sce, useReducedDim = "scanpyPCA", useDims = 1:2, nTop = 6, 
                    labelSize = 8, titleSize = 24, axisTextSize = 24)
dev.off()

pdf("scanpy_marker_dotPlot.pdf", width = 6, height = 5)
dp <- plotScanpyMarkerGenesDotPlot(sce, groupBy = "Scanpy_louvain_0.8", 
                                  nGenes = 4, vmin = -2, vmax = 2)
dp$show_colorbar <- FALSE
dp$show_size_legend <- FALSE
dp$largest_dot <- 120
dp$smallest_dot
dp$show()
dev.off()
  
# Downstream Analysis ####
sce <- runFindMarker(sce, useAssay = "logcounts", method = "wilcox",
                     cluster = "cluster")
sce <- runSingleR(sce, useAssay = "logcounts", level = "fine")
sce <- importGeneSetsFromMSigDB(sce, categoryIDs = "H", species = "Homo sapiens")
sce <- runVAM(sce, geneSetCollectionName = "H", useAssay = "logcounts")
hallmark <- "HALLMARK_INFLAMMATORY_RESPONSE"
sce <- runTSCAN(sce, useReducedDim = "PCA", cluster = "cluster")

## plot FindMarker Heatmap ####
pdf("findmarkerHeatmap.pdf", width = 3.6, height = 3)
plotFindMarkerHeatmap(sce, topN = 5, log2fcThreshold = 0.5,
                      fdrThreshold = 0.05, minClustExprPerc = 0.5,
                      maxCtrlExprPerc = 0.5, minMeanExpr = 0, rowLabel = FALSE,
                      row_title_gp = grid::gpar(fontsize = 12),
                      column_title_gp = grid::gpar(fontsize = 12),
                      show_heatmap_legend = FALSE)
dev.off()
## Plot singleR ####
pdf("singleR.pdf", width = 10, height = 9)
plotSCEDimReduceColData(sce, colorBy = "SingleR_hpca_fine_pruned.labels",
                        reducedDimName = "UMAP", baseSize = 28,
                        labelClusters = FALSE) +
    theme(legend.position = "none")
dev.off()
## Plot VAM ####
pdf("VAM.pdf", width = 10, height = 9)
plotPathway(sce, resultName = "VAM_H_CDF", geneset = hallmark,
            groupBy = "cluster", titleSize = 28, axisSize = 28,
            axisLabelSize = 28)
dev.off()
## Plot TSCAN ####
pdf("TSCAN.pdf", width = 10, height = 9)
plotTSCANResults(sce, useReducedDim = "UMAP") +
    theme(text = element_text(size = 28),
          axis.text = element_text(size = 28),
          legend.position = "none")
dev.off()

# Now that all vectorized sub-figures are generated, we use software "Inkscape"
# to put things together.


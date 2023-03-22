# IMPORTANT! RESTART R SESSION BEFORE EVERY RUN #

library(peakRAM)
library(dplyr)

if (!requireNamespace("singleCellTK"))
    devtools::install_github("compbiomed/singleCellTK@devel")
library(singleCellTK)
setwd("/restricted/projectnb/camplab/home/wangych/work/git/SCTK_benchmarking/")
# Only need to change this line for running different dataset
dataset <- "pbmc68k"

# Be sure to have raw data SCE object stored as flat files
#sce <- importExampleData(dataset)
#exportSCEtoFlatFile(sce, paste0(dataset, "_flat/"))

assayFile <- paste0(dataset, "_flat/assays/SCE_counts.mtx.gz")
annotFile <- paste0(dataset, "_flat/SCE_cellData.txt.gz")
featureFile <- paste0(dataset, "_flat/SCE_featureData.txt.gz")

steps <- factor(c("Import", "QC", "Filter", "Normalization","HVG",
                  "PCA", "UMAP", "cluster", "marker"), 
                levels = c("Import", "QC", "Filter", "Normalization",
                           "HVG", "PCA", "UMAP", "cluster", "marker"))

gc()

bm <- peakRAM(
    sce <- importFromFiles(assayFile = assayFile, 
                           annotFile = annotFile, featureFile = featureFile, 
                           annotFileHeader = TRUE, annotFileSep = ",", 
                           featureHeader = TRUE, featureSep = ","),
    sce <- runCellQC(sce, sample = NULL,
                     algorithms = c("QCMetrics"),
                     mitoRef = "human", mitoIDType = "symbol", 
                     mitoGeneLocation = "rownames", seed = 12345),
    sce <- subsetSCECols(sce, colData = c("total > 600", 
                                          "detected > 300",
                                          "mito_percent < 5")),
    sce <- runNormalization(sce, useAssay = "counts", outAssayName = "logcounts",
                            normalizationMethod = "logNormCounts"),
    {
        sce <- runModelGeneVar(sce, useAssay = "logcounts")
        sce <- setTopHVG(sce, method = "modelGeneVar", hvgNumber = 2000, featureSubsetName = "hvg2000")
    },
    sce <- scaterPCA(sce, useFeatureSubset = "hvg2000", seed = 12345),
    sce <- runUMAP(sce, useReducedDim = "PCA", initialDims = 10, seed = 12345),
    sce <- runScranSNN(sce, useReducedDim = "PCA", clusterName = "cluster", nComp = 10, seed = 12345),
    sce <- runFindMarker(sce, useAssay = "logcounts", method = "wilcox", cluster = "cluster")
)
bm <- bm %>% 
    mutate(Total_RAM_cumulative = cumsum(Total_RAM_Used_MiB),
           Step = steps, Dataset = dataset)
write.csv(bm, paste0(dataset, "_benchmark.csv"))



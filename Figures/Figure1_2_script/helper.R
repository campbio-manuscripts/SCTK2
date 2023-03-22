library(ggplot2)
library(singleCellTK)

plotScanpyPCLoading <- function(
    inSCE, 
    useReducedDim = "scanpyPCA", 
    useDims = seq(3),
    nTop = 10,
    nBottom = nTop,
    labelSize = 2,
    axisTextSize = 16,
    titleSize = 10) 
{
    pc <- reducedDim(inSCE, useReducedDim)
    loading <- attributes(p)[["rotation"]]
    if (is.null(loading)) 
        stop("No loading information found in the specified reducedDim")
    loading <- loading[,useDims]
    plotList <- list()
    for (j in seq(ncol(loading))) {
        loading.j <- data.frame(loading = sort(loading[,j], decreasing = TRUE))
        showIdx <- c(seq(nTop), 
                     seq(nrow(loading.j) - nBottom + 1, 
                         nrow(loading.j))
                     )
        loading.j <- loading.j[showIdx, , drop = FALSE]
        loading.j$name <- rownames(loading.j)
        loading.j$ranking <- seq(nrow(loading.j))
        insert <- nrow(loading.j) + 1
        loading.j[insert, "loading"] <- 
          (max(loading.j$loading) + min(loading.j$loading))/2
        loading.j[insert, "name"] <- "..."
        loading.j[insert, "ranking"] <- nTop + 0.5
        
        plotList[[j]] <- ggplot(loading.j, 
                                aes(x = ranking, y = loading, label = name)) +
            geom_text(angle = 90, size = labelSize) +
            ggtitle(paste0("PC", j)) + 
            ylim(min(loading.j$loading) - 0.035, max(loading.j$loading) + 0.035) +
            theme(panel.border = element_rect(fill = NA),
                  panel.background = element_blank(),
                  axis.text = element_text(size = axisTextSize - 2),
                  axis.title = element_text(size = axisTextSize),
                  title = element_text(size = titleSize),
                  axis.text.x = element_blank())
            
    }
    cowplot::plot_grid(plotlist = plotList, align = "hv", axis = "tblr")
}

#plotScanpyPCLoading(sce, useDims = 1:2, labelSize = 4, titleSize = 20, nTop = 5)


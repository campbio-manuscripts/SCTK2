library(ggplot2)
if (!require("ggbreak")) install.packages("ggbreak")
if (!require("ggforce")) install.packages("ggforce")
library(ggbreak)
library(ggforce)
# Only need to change this line for adding dataset
setwd("/Users/wangych/work/camplab/sctk/manuscript/Figure5")
source("facetZoom2.R")
dataOrder <- c("pbmc6k", "pbmc68k", "immune100k", "immune300k")


dfList <- lapply(dataOrder, function(x){
    read.csv(paste0(x, "_benchmark.csv"))
})
bm.df <- do.call(rbind, dfList)
# Try to shorten the word so maybe I can enlarge the text in the figure
bm.df$Step[bm.df$Step == "Import"] <- "import"
bm.df$Step[bm.df$Step == "Filter"] <- "filter"
bm.df$Step[bm.df$Step == "Normalization"] <- "norm"
stepOrder <- c("import", "QC", "filter", "norm", "HVG", 
               "PCA", "UMAP", "cluster", "marker")
bm.df$Step <- factor(bm.df$Step, levels = stepOrder)
bm.df$Dataset <- factor(bm.df$Dataset, levels = dataOrder)

trim <- function(x, min = NULL, max = NULL) {
    if (!is.null(min)) x <- ifelse(x < min, min, x)
    if (!is.null(max)) x <- ifelse(x > max, max, x)
    return(x)
}
axis.text.fontsize <- 13


p1 <- ggplot() +
    geom_bar(aes(x = factor(Step, levels=levels(Step)), 
                 y = Peak_RAM_Used_MiB/1000, 
                 fill = Dataset), 
             data = bm.df, position="dodge", stat="identity") +
    xlab("Steps") + ylab("Peak RAM Used (GB)") +
    scale_fill_discrete("Dataset") +
    theme_classic() + 
    theme(axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = 70, vjust = 0.9, hjust=1, size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 13))

p2 <- ggplot() +
    geom_line(aes(x = factor(Step, levels=levels(Step)), 
                  y = Total_RAM_cumulative/1000, 
                  color = Dataset, 
                  group = Dataset), 
              data = bm.df) +
    xlab("Steps") + ylab("RAM of object (GB)") +
    scale_fill_discrete("Dataset") +
    theme_classic() +  
    theme(axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = 70, vjust = 0.9, hjust=1, size = 12))

p3 <- ggplot() +
    geom_bar(aes(x = factor(Step, levels=levels(Step)), 
                 y = Elapsed_Time_sec/60,
                 fill = Dataset), 
             data = bm.df, position="dodge", stat="identity") +
    xlab("Steps") + ylab("Time Elapsed (min)") +
    scale_fill_discrete("Dataset") +
    theme_classic() +  
    theme(axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = 70, vjust = 0.9, hjust=1, size = 12),
          strip.background = element_rect(fill = "grey85", 
                                          colour = NA), ) +
    facet_zoom2(ylim = c(0, 7), zoom.size = 2)

legend <- cowplot::get_legend(p1)
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
up <- cowplot::plot_grid(p2, p1,
                         labels = c("A", "B"),
                         ncol = 2)
down <- cowplot::plot_grid(p3, legend, 
                   labels = c("C", ""),
                   ncol = 2, rel_widths = c(4, 1))
cowplot::plot_grid(up,
                   down,
                   ncol = 1)


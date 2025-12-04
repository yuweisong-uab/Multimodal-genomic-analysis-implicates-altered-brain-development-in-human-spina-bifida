set.seed(1234)

# Core Single-Cell Analysis & Visualization
library(Seurat)
library(Signac)
library(SeuratData)
library(SeuratWrappers)
library(monocle3)
library(Nebulosa)
library(clustree)

# Genomic Data and Annotation
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(biovizBase)
library(ensembldb)
library(biomaRt)

# Data Processing & Statistical Analysis
library(Matrix)
library(metap)
library(multtest)
library(glmGamPoi)
library(future)
library(vegan)

# Visualization Tools
library(ggplot2)
library(ggpubr)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggrepel)
library(scales)
library(dplyr)

# File Handling & Web Access
library(RCurl)
library(openxlsx)

# Functional Annotation & Enrichment Analysis
library(AnnotationHub)
library(org.Hs.eg.db)
library(clusterProfiler)

# Single-Cell Proteomics & Others
library(SCP)
library(genekitr)

# Font Support
library(showtext)
library(extrafont)
loadfonts(device="pdf")

# Python Interface
library(reticulate)

# Optional raster setting
raster = FALSE

#Functions

style_umap_axes <- function(p, obj, scale_len = 0.2, axis_size = 5,
                            legend_key = 0.3, legend_text = 7, legend_title = 8) {
  
  umap <- Embeddings(obj, "umap")
  x_range <- range(umap[,1])
  y_range <- range(umap[,2])
  
  p +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = c(0.05, 0.95),
      legend.justification = c(0, 1),
      legend.background = element_blank(),
      legend.key.size = unit(legend_key, "cm"), 
      legend.text = element_text(size = legend_text),
      legend.title = element_text(size = legend_title)
    ) +
    geom_segment(
      aes(x = x_range[1],
          xend = x_range[1] + diff(x_range)*scale_len,
          y = y_range[1],
          yend = y_range[1]),
      arrow = arrow(type="closed", length = unit(0.22,"cm")),
      linewidth = 0.6
    ) +
    geom_segment(
      aes(x = x_range[1],
          xend = x_range[1],
          y = y_range[1],
          yend = y_range[1] + diff(y_range)*scale_len),
      arrow = arrow(type="closed", length = unit(0.22,"cm")),
      linewidth = 0.6
    ) +
    annotate("text",
             x = x_range[1] + diff(x_range)*scale_len/2,
             y = y_range[1] - diff(y_range)*0.05,
             label = "UMAP1",
             size = axis_size) +
    annotate("text",
             x = x_range[1] - diff(x_range)*0.05,
             y = y_range[1] + diff(y_range)*scale_len/2,
             label = "UMAP2",
             angle = 90,
             size = axis_size)
}

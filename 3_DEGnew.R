setwd("~/Library/CloudStorage/OneDrive-UAB-TheUniversityofAlabamaatBirmingham/1_ChongLab/Project4_Andy_hydrocephalus/DuyPhan_genecheck/")

source(file = "scripts/Plot_colorPaletters.R")
source(file = "scripts/load_libraries.R")
plan("multisession", workers = 4) 
options(future.globals.maxSize = 50000 * 1024^2) # set 50G RAM

load("output/RData/3_integrate_HCPUAB.RData")
my_cols <- c('Astrocytes'='#8F7700','OPC'='#3B3B3B','Oligodendrocytes'='#EFC001','Microglia'='#CD534C','Excitatory neurons'='#59753E','Inhibitory neurons'='#394E2F')
Idents(HCPUAB)<-HCPUAB$celltype
table(HCPUAB$celltype,HCPUAB$sampleID)
cluster.all.markers <- FindAllMarkers(HCPUAB, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25,verbose = T,assay = 'RNA')
View(cluster.all.markers)
write.csv(cluster.all.markers,file="output/summary/all_cell_markers.csv")
cluster.sig.markers <- cluster.all.markers[which(cluster.all.markers$p_val_adj<0.1),]
write.csv(cluster.sig.markers, file = "output/summary/HCPUAB_cluster.sig.markers.csv")

top.genes <- cluster.sig.markers %>%group_by(cluster) %>%slice_max(order_by = avg_log2FC, n = 15)
write.csv(top.genes, file = "output/summary/HCPUAB_cluster.sig.markers.top10.csv")

#GeneList_tocheck <- read.xlsx("pnas.2106844118.sd01.xlsx")
#GeneList_tocheck <- GeneList_tocheck$Symbol
# 
# ##DEG in general
# HCPUAB <- RunDEtest(srt = HCPUAB, group_by = "celltype",  only.pos = F,fc.threshold = 1)
# VolcanoPlot(srt = HCPUAB, group_by = "celltype")
# DEGs <- HCPUAB@tools$DEtest_celltype$AllMarkers_wilcox
# DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]
# #%>% group_by(group1) %>% top_n(n = 50, wt = avg_log2FC)
# # Annotate features with transcription factors and surface proteins
# DEGs$group1<-factor(DEGs$group1,levels=c("Astrocytes",  "OPC", "Oligodendrocytes","Inhibitory neurons", 'Excitatory neurons',"Microglia"))
# cluster.sig.markers$cluster <-factor(cluster.sig.markers$cluster,levels=c("Astrocytes",  "OPC", "Oligodendrocytes","Inhibitory neurons", 'Excitatory neurons',"Microglia"))
# 
# HCPUAB <- AnnotateFeatures(HCPUAB, species = "Homo_sapiens", db = c("TF", "CSPA"))
# ht <- FeatureHeatmap(
#   srt = HCPUAB, group.by = "celltype", features = cluster.sig.markers$gene, 
#   species = "Homo_sapiens", #feature_annotation = c("TF", "CSPA"),
#   #anno_terms = TRUE,
#   #db = c("GO_BP", "KEGG"), 
#   #feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
#   height = 5, width = 4
# )
# ht$plot
# 
# 
# 
# 
# 
# 
# 
# ht <- GroupHeatmap(
#   srt = HCPUAB,
#   features = top.genes$gene,
#   group.by="celltype",
#   group_palcolor = my_cols,
#   #show_row_names = TRUE, row_names_side = "left",
#   add_dot = TRUE, add_reticle = TRUE
# )
# 
# print(ht$plot)library(Seurat)
library(Seurat)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

## ========= 0. Settings =========
pair_name <- "Pair_HCPUAB014"
ident.1   <- "HCPUAB014"
ident.2   <- "control_y10"
ct        <- "OPC"     # one celltype to test

deg_dir <- "output/summary/DEG/"

non_coding_keywords <- c(
  "^AC\\d+", "^AL\\d+", "-AS\\d*$", "^LINC", "^MIR", "^SNO", "^sno",
  "^SNR", "^RPS", "^RPL", "^ENSG", "-DT\\d*$"
)

significance_threshold <- 0.05
contrast_label <- paste(ident.1, "vs", ident.2)

## ========= 1. Subset to this pair =========
obj_pair <- subset(HCPUAB, subset = orig.ident %in% c(ident.1, ident.2))
obj_pair$orig.ident <- factor(obj_pair$orig.ident, levels = c(ident.1, ident.2))

## check cell numbers per celltype
cell_counts <- table(obj_pair$orig.ident, obj_pair$celltype)
print(cell_counts)

if (!ct %in% colnames(cell_counts)) {
  stop("Cell type ", ct, " not present in this pair.")
}

min_ct <- min(cell_counts[, ct])
if (min_ct < 20) {
  stop("Cell type ", ct, " has too few cells (<20) in at least one group: ", min_ct)
}

## ========= 2. Subset to this celltype =========
srt_ct <- subset(obj_pair, subset = celltype == ct)

if (length(unique(srt_ct$orig.ident)) < 2) {
  stop("Only one orig.ident group present in ", ct, " for this pair.")
}

## ========= 3. RunDEtest for volcano + heatmap =========
srt_ct <- RunDEtest(
  srt = srt_ct,
  group_by = "orig.ident",
  fc.threshold = 1,
  only.pos = FALSE
)

# Volcano plot
p_vol <- VolcanoPlot(srt = srt_ct, group_by = "orig.ident")
ggsave(
  filename = file.path(
    deg_dir,
    paste0("Volcano_", pair_name, "_", gsub(" ", "_", ct), "_Volcano.pdf")
  ),
  plot = p_vol, width = 8, height = 5
)

## ---- 3.1 Heatmap based on DEtest result ----
DEGs_DEtest <- srt_ct@tools$DEtest_orig.ident$AllMarkers_wilcox
DEGs_DEtest <- DEGs_DEtest[with(DEGs_DEtest, avg_log2FC > 1 & p_val_adj < 0.05), ]

DEGs_filtered <- DEGs_DEtest %>%
  filter(!str_detect(gene, paste(non_coding_keywords, collapse = "|")))

if (nrow(DEGs_filtered) > 0) {
  top.30genes <- DEGs_filtered %>%
    group_by(group1) %>%
    dplyr::top_n(n = 30, wt = avg_log2FC)
  
  ht <- FeatureHeatmap(
    srt = srt_ct,
    features = top.30genes$gene,
    feature_split = top.30genes$group1,
    group.by = "orig.ident",
    species = "Homo_sapiens",
    row_names_side = "right",
    show_row_names = TRUE,
    column_title = paste0(contrast_label, " — ", ct),
    height = 10,
    width  = 5
  )
  
  pdf(
    file.path(
      deg_dir,
      paste0(pair_name, "_", gsub(" ", "_", ct), "_DEG_heatmap.pdf")
    ),
    height = 13, width = 12
  )
  print(ht$plot)
  dev.off()
} else {
  message("No filtered DEGs (for heatmap) in ", ct)
}

## ========= 4. FindMarkers for DEG / GO / KEGG =========
Idents(srt_ct) <- srt_ct$orig.ident
DE_all <- FindMarkers(
  srt_ct,
  ident.1 = ident.1,
  ident.2 = ident.2,
  verbose = TRUE,
  assay   = "RNA"
)

sig_dge.all  <- subset(DE_all, p_val_adj < 0.05)
if (nrow(sig_dge.all) == 0) {
  stop("No significant DEGs (FindMarkers) in ", ct, " — cannot run GO/KEGG.")
}

sig_dge.up   <- subset(sig_dge.all, avg_log2FC >  0.15)
sig_dge.down <- subset(sig_dge.all, avg_log2FC < -0.15)

sig_dge.up   <- sig_dge.up[order(sig_dge.up$avg_log2FC, decreasing = TRUE), ]
sig_dge.down <- sig_dge.down[order(sig_dge.down$avg_log2FC, decreasing = FALSE), ]

base_name <- paste0(pair_name, "_", gsub(" ", "_", ct))

write.csv(sig_dge.all,
          file.path(deg_dir, paste0(base_name, "_DEG_all_sig.csv")))
write.csv(sig_dge.up,
          file.path(deg_dir, paste0(base_name, "_DEG_sig_up.csv")))
write.csv(sig_dge.down,
          file.path(deg_dir, paste0(base_name, "_DEG_sig_down.csv")))

## ========= 5. GO BP enrichment (up / down) =========
GO_BP_up_df   <- NULL
GO_BP_down_df <- NULL

# up
if (nrow(sig_dge.up) > 0) {
  GO_BP_up <- enrichGO(
    gene          = rownames(sig_dge.up),
    OrgDb         = "org.Hs.eg.db",
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05
  )
  if (!is.null(GO_BP_up@result) && nrow(GO_BP_up@result) > 0) {
    GO_BP_up_df <- as.data.frame(GO_BP_up@result)
    GO_BP_up_df$regulation <- "genes increased in patient"
    GO_BP_up_df <- GO_BP_up_df[order(GO_BP_up_df$pvalue), ]
    GO_BP_up_df_sig <- GO_BP_up_df %>%
      filter(pvalue < significance_threshold)
    write.csv(
      GO_BP_up_df_sig,
      file.path(deg_dir, paste0(base_name, "_GO_BP_up_sig.csv"))
    )
  }
}

# down
if (nrow(sig_dge.down) > 0) {
  GO_BP_down <- enrichGO(
    gene          = rownames(sig_dge.down),
    OrgDb         = "org.Hs.eg.db",
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05
  )
  if (!is.null(GO_BP_down@result) && nrow(GO_BP_down@result) > 0) {
    GO_BP_down_df <- as.data.frame(GO_BP_down@result)
    GO_BP_down_df$regulation <- "genes decreased in patient"
    GO_BP_down_df <- GO_BP_down_df[order(GO_BP_down_df$pvalue), ]
    GO_BP_down_df_sig <- GO_BP_down_df %>%
      filter(pvalue < significance_threshold)
    write.csv(
      GO_BP_down_df_sig,
      file.path(deg_dir, paste0(base_name, "_GO_BP_down_sig.csv"))
    )
  }
}

# combined GO_BP barplot
if (!is.null(GO_BP_up_df) | !is.null(GO_BP_down_df)) {
  combined_GO_BP <- dplyr::bind_rows(
    if (!is.null(GO_BP_up_df))   head(GO_BP_up_df,   10),
    if (!is.null(GO_BP_down_df)) head(GO_BP_down_df, 10)
  )
  
  if (nrow(combined_GO_BP) > 0) {
    combined_GO_BP <- combined_GO_BP %>%
      mutate(
        log_pvalue = -log10(pvalue),
        log_pvalue = ifelse(
          regulation == "genes decreased in patient",
          -log_pvalue, log_pvalue
        )
      ) %>%
      arrange(regulation, log_pvalue) %>%
      mutate(Description = factor(Description, levels = unique(Description)))
    
    write.csv(
      combined_GO_BP,
      file.path(deg_dir, paste0(base_name, "_GO_BP_combined.csv")),
      row.names = FALSE
    )
    
    p_go <- ggplot(
      combined_GO_BP,
      aes(x = reorder(Description, log_pvalue),
          y = log_pvalue, fill = regulation)
    ) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(
        values = c(
          "genes increased in patient" = "#D2413C",
          "genes decreased in patient" = "#367EB9"
        ),
        name = NULL
      ) +
      scale_y_continuous(labels = function(x) abs(x)) +
      labs(
        x = "",
        y = "-log10(p-value)",
        title = paste0(
          "GO enrichment of differentially expressed genes in ", ct, "\n",
          ident.1, " patient versus ", ident.2
        )
      ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8),
        axis.title  = element_text(color = "black", size = 10),
        plot.title  = element_text(hjust = 0.5, size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(fill = NA, color = "black", size = 1),
        plot.background  = element_rect(fill = "white", color = NA),
        axis.line        = element_line(color = "black"),
      )
    
    ggsave(
      file.path(deg_dir, paste0("legend.pdf")),
      p_go, width = 8, height = 4
    )
  }
}

## ========= 6. KEGG enrichment (up / down) =========
KEGG_up_df   <- NULL
KEGG_down_df <- NULL

# up
if (nrow(sig_dge.up) > 0) {
  genelist_up <- clusterProfiler::bitr(
    rownames(sig_dge.up),
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = "org.Hs.eg.db"
  )
  if (!is.null(genelist_up) && nrow(genelist_up) > 0) {
    KEGG_up <- clusterProfiler::enrichKEGG(
      gene     = genelist_up$ENTREZID,
      organism = "hsa"
    )
    if (!is.null(KEGG_up@result) && nrow(KEGG_up@result) > 0) {
      KEGG_up_df <- as.data.frame(KEGG_up@result)
      KEGG_up_df$regulation <- "genes increased in patient"
      KEGG_up_df <- KEGG_up_df[order(KEGG_up_df$pvalue), ]
      KEGG_up_df_sig <- KEGG_up_df %>%
        dplyr::filter(pvalue < significance_threshold)
      write.csv(
        KEGG_up_df_sig,
        file.path(deg_dir, paste0(base_name, "_KEGG_up_sig.csv"))
      )
    }
  }
}

# down
if (nrow(sig_dge.down) > 0) {
  genelist_down <- clusterProfiler::bitr(
    rownames(sig_dge.down),
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = "org.Hs.eg.db"
  )
  if (!is.null(genelist_down) && nrow(genelist_down) > 0) {
    KEGG_down <- clusterProfiler::enrichKEGG(
      gene     = genelist_down$ENTREZID,
      organism = "hsa"
    )
    if (!is.null(KEGG_down@result) && nrow(KEGG_down@result) > 0) {
      KEGG_down_df <- as.data.frame(KEGG_down@result)
      KEGG_down_df$regulation <- "genes decreased in patient"
      KEGG_down_df <- KEGG_down_df[order(KEGG_down_df$pvalue), ]
      KEGG_down_df_sig <- KEGG_down_df %>%
        dplyr::filter(pvalue < significance_threshold)
      write.csv(
        KEGG_down_df_sig,
        file.path(deg_dir, paste0(base_name, "_KEGG_down_sig.csv"))
      )
    }
  }
}

# combined KEGG barplot
if (!is.null(KEGG_up_df) | !is.null(KEGG_down_df)) {
  combined_KEGG <- dplyr::bind_rows(
    if (!is.null(KEGG_up_df))   head(KEGG_up_df,   10),
    if (!is.null(KEGG_down_df)) head(KEGG_down_df, 10)
  )
  
  if (nrow(combined_KEGG) > 0) {
    combined_KEGG <- combined_KEGG %>%
      dplyr::mutate(
        log_pvalue = -log10(pvalue),
        log_pvalue = ifelse(
          regulation == "genes decreased in patient",
          -log_pvalue, log_pvalue
        )
      ) %>%
      dplyr::arrange(regulation, log_pvalue) %>%
      dplyr::mutate(Description = factor(Description, levels = unique(Description)))
    
    write.csv(
      combined_KEGG,
      file.path(deg_dir, paste0(base_name, "_KEGG_combined.csv")),
      row.names = FALSE
    )
    
    p_kegg <- ggplot(
      combined_KEGG,
      aes(x = reorder(Description, log_pvalue),
          y = log_pvalue, fill = regulation)
    ) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(
        values = c(
          "genes increased in patient" = "#D2413C",
          "genes decreased in patient" = "#367EB9"
        ),
        name = NULL
      ) +
      scale_y_continuous(labels = function(x) abs(x)) +
      labs(
        x = "",
        y = "-log10(p-value)",
        title = paste0(
          "KEGG enrichment of differentially expressed genes in ", ct, "\n",
          ident.1, " patient versus ", ident.2
        )
      ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8),
        axis.title  = element_text(color = "black", size = 10),
        plot.title  = element_text(hjust = 0.5, size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(fill = NA, color = "black", size = 1),
        plot.background  = element_rect(fill = "white", color = NA),
        axis.line        = element_line(color = "black"),
        legend.position  = "none"
      )
    
    ggsave(
      file.path(deg_dir, paste0(base_name, "_KEGG_barplot.pdf")),
      p_kegg, width = 8, height = 4
    )
  }
}




setwd("~/Library/CloudStorage/OneDrive-UAB-TheUniversityofAlabamaatBirmingham/1_ChongLab/Project4_Andy_hydrocephalus/DuyPhan_genecheck/")

source(file = "scripts/Plot_colorPaletters.R")
source(file = "scripts/load_libraries.R")

load("output/RData/3_integrate_HCPUAB.RData")
# HCPUAB[["RNA3"]] <- as(object = HCPUAB[["RNA"]], Class = "Assay")
# DefaultAssay(HCPUAB) <- "RNA3"
# HCPUAB[["RNA"]] <- NULL
# HCPUAB[["ATAC"]] <- NULL
#HCPUAB <- RenameAssays(object = HCPUAB, RNA3 = 'RNA')
my_cols <- c('Astrocytes'='#8F7700','OPC'='#3B3B3B','Oligodendrocytes'='#EFC001','Microglia'='#CD534C','Excitatory neurons'='#59753E','Inhibitory neurons'='#394E2F')

#check cell number
cellnumber_sample_Etiology<-table(HCPUAB@meta.data$orig.ident)
write.csv(cellnumber_sample_Etiology,file="output/summary/cellnumber_sample_Etiology.csv")

#set factor
HCPUAB$SampleID<-factor(HCPUAB$SampleID,levels=c("1M control","1M patient","10Y control","10Y patient"))
HCPUAB$Condition <-factor(HCPUAB$Condition,levels=c("Control","Patient"))
HCPUAB$celltype <- factor(HCPUAB$celltype, 
                                   levels = c("Astrocytes",  "OPC", "Oligodendrocytes","Inhibitory neurons", 'Excitatory neurons',"Microglia"))


# #### find all markers
cluster.all.markers <- FindAllMarkers(HCPUAB, only.pos = TRUE, min.pct = 0.1,group.by = c("seurat_clusters"),logfc.threshold = 0.25,verbose = T,assay = 'RNA')
view(cluster.all.markers)
write.csv(cluster.all.markers,file="output/summary/all_cell_markers.csv")
cluster.sig.markers <- cluster.all.markers[which(cluster.all.markers$p_val_adj<0.1),]
write.csv(cluster.sig.markers, file = "output/summary/HCPUAB_cluster.sig.markers.csv")
top.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
write.csv(top.genes, file = "output/summary/HCPUAB_cluster.sig.markers.top10.csv") 

# #use GPTcellType for annotation
Sys.setenv(OPENAI_API_KEY = '')
library(GPTCelltype)
library(openai)
res <- gptcelltype(cluster.sig.markers, tissuename = 'human brain',model = 'gpt-4')

cat(res)



# Visualize cell type annotation on UMAP
Idents(HCPUAB)<-HCPUAB$seurat_clusters
DimPlot(HCPUAB,label = T)
table(HCPUAB@meta.data$seurat_clusters)
cluster.label <- HCPUAB@meta.data$seurat_clusters
cluster.label <- gsub("^0$", "Astrocytes", cluster.label)
cluster.label <- gsub("^1$", "OPC", cluster.label)
cluster.label <- gsub("^2$", "Oligodendrocytes", cluster.label)
cluster.label <- gsub("^3$", "Microglia", cluster.label)
cluster.label <- gsub("^4$", "Astrocytes", cluster.label)
cluster.label <- gsub("^5$", "Astrocytes", cluster.label)
cluster.label <- gsub("^6$", 'Inhibitory neurons', cluster.label)
cluster.label <- gsub("^7$", 'Excitatory neurons', cluster.label) 
cluster.label <- gsub("^8$", 'Excitatory neurons', cluster.label)
cluster.label <- gsub("^9$", 'Inhibitory neurons', cluster.label)
cluster.label <- gsub("^10$", "Oligodendrocytes", cluster.label)
cluster.label <- gsub("^11$", 'Excitatory neurons', cluster.label) 
cluster.label <- gsub("^12$", "Astrocytes", cluster.label)
cluster.label <- gsub("^13$", 'Inhibitory neurons', cluster.label)
cluster.label <- gsub("^14$", "Astrocytes", cluster.label)
cluster.label <- gsub("^15$", "Oligodendrocytes", cluster.label)
cluster.label <- gsub("^16$", 'Excitatory neurons', cluster.label)
cluster.label <- gsub("^17$", "OPC", cluster.label)
cluster.label <- gsub("^18$", "Oligodendrocytes", cluster.label)

table(HCPUAB$seurat_clusters)
HCPUAB<- AddMetaData(HCPUAB, cluster.label, col.name = "celltype")
DimPlot(HCPUAB,group.by='celltype',label = TRUE)
#HCPUAB<- subset(HCPUAB, subset = seurat_clusters %in% c(0:18))
table(HCPUAB$celltype)
DimPlot(HCPUAB,group.by='seurat_clusters',label = TRUE)
DimPlot(HCPUAB,group.by='celltype',label = TRUE)
HCPUAB@meta.data$celltype<- factor(HCPUAB@meta.data$celltype, 
                                   levels = c("Astrocytes",  "OPC", "Oligodendrocytes","Inhibitory neurons", 'Excitatory neurons',"Microglia"))

#my_cols <- c('Astrocytes'='#E7C359','OPC'='#93ADC0','Oligodendrocytes'='#2C426F','Microglia'='#AA5F48','Excitatory neurons'='#59753E','Inhibitory neurons'='#394E2F')
my_cols <- c('Astrocytes'='#8F7700','OPC'='#3B3B3B','Oligodendrocytes'='#EFC001','Microglia'='#CD534C','Excitatory neurons'='#59753E','Inhibitory neurons'='#394E2F')


DimPlot(HCPUAB,group.by='celltype',label = TRUE,cols = my_cols)

Idents(HCPUAB)<-HCPUAB$celltype


##提取umap数据
HCPUAB[['celltype']] = HCPUAB@active.ident
umap = HCPUAB@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(celltype = HCPUAB@meta.data$celltype)

##计算标签中心位置区域
celltypepos <- umap %>%
  group_by(celltype) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2))


##integrated umap in general 
#pdf("output/summary/HCPUAB_annotation.pdf", height = 12, width=12,family = "Arial")
p <-  DimPlot(HCPUAB, reduction = "umap", label = FALSE, cols = my_cols,group.by = "celltype") + 
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 12,hjust = 0.5),  # Center the x-axis label
        axis.title.y = element_text(size = 12,vjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  NoLegend() +
  geom_label_repel(aes(x = umap_1,y = umap_2,label = celltype,color = celltype), 
                   fontface = "bold",
                   data = celltypepos,
                   box.padding = 0.5,
                   size = 3.5)+
  labs(x = "UMAP1", y = "UMAP2")+
  ggtitle("Integrated dataset")
p
ggsave("output/summary/HCPUAB_annotation.pdf", height = 5, width=5,family = "Arial", p, device = cairo_pdf)
#dev.off()

## UMAP condition no label
p1 <- DimPlot(HCPUAB, reduction = "umap", label = FALSE,cols = my_cols,
              group.by = "celltype", split.by = "Condition", ncol = 2) +
  NoLegend()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(hjust = 0.5,size = 12),
        axis.title.y = element_text(vjust = 0.5,size = 12),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 12)) +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("")

ggsave("output/summary/HCPUAB_annotation_condition_nolabel.pdf",
       p1, width = 12, height = 8, device = cairo_pdf, family = "Arial")


## UMAP condition with label
p2 <- DimPlot(HCPUAB, reduction = "umap", label = FALSE,cols = my_cols,
              group.by = "celltype", split.by = "Condition", ncol = 2) +
  NoLegend()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(hjust = 0.5,size = 12),
        axis.title.y = element_text(vjust = 0.5,size = 12),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 12)) +
  geom_label_repel(aes(x = umap_1,y = umap_2,label = celltype,color = celltype),
                   data = celltypepos, fontface = "bold",
                   box.padding = 0.5, size = 5) +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("")

ggsave("output/summary/HCPUAB_annotation_condition_label.pdf",
       p2, width = 12, height = 8, device = cairo_pdf, family = "Arial")


## UMAP sample no label
p3 <- DimPlot(HCPUAB, reduction = "umap", label = FALSE,cols = my_cols,
              group.by = "celltype", split.by = "SampleID", ncol = 4) +
  NoLegend()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(hjust = 0.5,size = 12),
        axis.title.y = element_text(vjust = 0.5,size = 12),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 12)) +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("")

ggsave("output/summary/HCPUAB_annotation_sample_nolabel.pdf",
       p3, width = 24, height = 8, device = cairo_pdf, family = "Arial")


## UMAP sample with label
p4 <- DimPlot(HCPUAB, reduction = "umap", label = FALSE,cols = my_cols,
              group.by = "celltype", split.by = "SampleID", ncol = 4) +
  NoLegend()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(hjust = 0.5,size = 12),
        axis.title.y = element_text(vjust = 0.5,size = 12),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 12)) +
  geom_label_repel(aes(x = umap_1,y = umap_2,label = celltype,color = celltype),
                   data = celltypepos, fontface = "bold",
                   box.padding = 0.5, size = 5) +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("")

ggsave("output/summary/HCPUAB_annotation_sample_label.pdf",
       p4, width = 24, height = 8, device = cairo_pdf, family = "Arial")


## Barplot sample vertical
p5 <- ggplot(HCPUAB@meta.data %>%
               mutate(SampleID = factor(SampleID,
                                        levels = rev(c("1M control","1M patient",
                                                       "10Y control","10Y patient"))))) +
  geom_bar(aes(y = SampleID, fill = celltype), position = position_fill()) +
  scale_fill_manual(values = my_cols) +
  scale_x_continuous(labels = scales::percent) +
  guides(fill = guide_legend(title = NULL))+
  theme_bw() +
  theme(axis.text.y  = element_text(colour = "black"),
        axis.text.x  = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(1.5,"line"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10)) +
  labs(y = "", x = "") +
  ggtitle("")

ggsave("output/summary/HCPUABcell_number_barplot_sample_vertical.pdf",
       p5, width = 7, height = 3, device = cairo_pdf, family = "Arial")


## Barplot sample horizontal
p6 <- ggplot(HCPUAB@meta.data %>%
               mutate(SampleID = factor(SampleID,
                                        levels = c("1M control","1M patient",
                                                   "10Y control","10Y patient")))) +
  geom_bar(aes(x = SampleID, fill = celltype), position = position_fill()) +
  scale_fill_manual(values = my_cols) +
  scale_y_continuous(labels = scales::percent) +
  guides(fill = guide_legend(title = NULL)) +
  theme_bw() +
  theme(axis.text.y  = element_text(colour = "black",size=12),
        axis.text.x  = element_text(colour = "black",size=12,
                                    angle = 30, hjust = 0.5,vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(2.5, "line"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10)) +
  labs(y = "", x = "") +
  ggtitle("")

ggsave("output/summary/HCPUABcell_number_barplot_sample_horizontal.pdf",
       p6, width = 7, height = 7, device = cairo_pdf, family = "Arial")


## Barplot condition vertical
p7 <- ggplot(HCPUAB@meta.data %>%
               mutate(Condition = factor(Condition,levels = c("Patient","Control")))) +
  geom_bar(aes(y = Condition, fill = celltype), position = position_fill()) +
  scale_fill_manual(values = my_cols) +
  scale_x_continuous(labels = scales::percent) +
  guides(fill = guide_legend(title = NULL))+
  theme_bw() +
  theme(axis.text.y  = element_text(colour = "black"),
        axis.text.x  = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(1.5,"line"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10)) +
  labs(y = "", x = "") +
  ggtitle("")

ggsave("output/summary/HCPUABcell_number_barplot_conditions_vertical.pdf",
       p7, width = 7, height = 2, device = cairo_pdf, family = "Arial")


## Barplot condition horizontal
p8 <- ggplot(HCPUAB@meta.data %>%
               mutate(Condition = factor(Condition,levels = c("Control","Patient")))) +
  geom_bar(aes(x = Condition, fill = celltype), position = position_fill()) +
  scale_fill_manual(values = my_cols) +
  scale_y_continuous(labels = scales::percent) +
  guides(fill = guide_legend(title = NULL)) +
  theme_bw() +
  theme(axis.text.y  = element_text( colour = "black",size=12),
        axis.text.x  = element_text( colour = "black",size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(2, "line"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10)) +
  labs(y = "", x = "") +
  ggtitle("")

ggsave("output/summary/HCPUABcell_number_barplot_conditions_horizontal.pdf",
       p8, width = 4, height = 5, device = cairo_pdf, family = "Arial")

#some published marker check
# #### find all markers
cluster.all.markers <- FindAllMarkers(HCPUAB, only.pos = TRUE, min.pct = 0.1,group.by = c("celltype"),logfc.threshold = 0.25,verbose = T,assay = 'RNA')
view(cluster.all.markers)
write.csv(cluster.all.markers,file="output/summary/all_celltype_markers.csv")
cluster.sig.markers <- cluster.all.markers[which(cluster.all.markers$p_val_adj<0.1),]
write.csv(cluster.sig.markers, file = "output/summary/HCPUAB_celltype.sig.markers.csv")
top.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
write.csv(top.genes, file = "output/summary/HCPUAB_celltype.sig.markers.top10.csv") 


# library(RColorBrewer)
# p1<-FeaturePlot(HCPUAB,features ="AQP4",cols = c("#D3D3D3","#8F7700"))+
#   labs(title = expression(Astrocytes~(italic("AQP4"))),x = NULL, y = NULL)
# p1 <- style_umap_axes(p1, HCPUAB)
# p2<-FeaturePlot(HCPUAB,features ="PDGFRA",cols = c("#D3D3D3","#3B3B3B")) +
#   labs(title = expression(OPCs~(italic("PDGFRA"))),x = NULL, y = NULL)
# p2 <- style_umap_axes(p2, HCPUAB)
# p3<-FeaturePlot(HCPUAB,features ="PLP1",cols = c("#D3D3D3","#EFC001"))+
#   labs(title = expression(Oligodendrocytes~(italic("PLP1"))),x = NULL, y = NULL)
# p3 <- style_umap_axes(p3, HCPUAB)
# 
# p4<-FeaturePlot(HCPUAB,features ="ADARB2",cols = c("#D3D3D3","#394E2F")) +
#   labs(title = expression(IN~(italic("ADARB2"))),x = NULL, y = NULL)
# p4 <- style_umap_axes(p4, HCPUAB)
# p5<-FeaturePlot(HCPUAB,features ="CCBE1",cols = c("#D3D3D3","#59753E")) +
#   labs(title = expression(EN~(italic("CCBE1"))),x = NULL, y = NULL)
# p5 <- style_umap_axes(p5, HCPUAB)
# p6<-FeaturePlot(HCPUAB,features ="PTPRC",cols = c("#D3D3D3","#CD534C")) +
#   labs(title = expression(Microglia~(italic("PTPRC"))),x = NULL, y = NULL)
# p6 <- style_umap_axes(p6, HCPUAB)
# pdf("output/summary/HCPUAB_celltype_FeaturePlot.pdf", height = 6, width =10 )
# (p1|p2|p3)/(p4|p5|p6)
# dev.off()

p1<-FeaturePlot(HCPUAB,features ="AQP4",cols = c("#D3D3D3","#8F7700"))+
  labs(title = expression(Astrocytes~(italic("AQP4"))),x = "UMAP1", y = "UMAP2")+
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 7))
p2<-FeaturePlot(HCPUAB,features ="PDGFRA",cols = c("#D3D3D3","#3B3B3B"))+
  labs(title = expression(OPCs~(italic("PDGFRA"))),x = "UMAP1", y = "UMAP2")+
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = c(0, 1),
      legend.justification = c(0, 1),
      legend.background = element_blank(),
      legend.key.size = unit(0.25, "cm"),
      legend.text = element_text(size = 7))

p3<-FeaturePlot(HCPUAB,features ="PLP1",cols = c("#D3D3D3","#EFC001"))+
  labs(title = expression(Oligodendrocytes~(italic("PLP1"))),x = "UMAP1", y = "UMAP2")+
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 7))
#03BF7D
p4<-FeaturePlot(HCPUAB,features ="ADARB2",cols = c("#D3D3D3","#394E2F")) +
  labs(title = expression(IN~(italic("ADARB2"))),x = "UMAP1", y = "UMAP2")+
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 7))
p5<-FeaturePlot(HCPUAB,features ="CCBE1",cols = c("#D3D3D3","#59753E")) +
  labs(title = expression(EN~(italic("CCBE1"))),x = "UMAP1", y = "UMAP2")+
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 7))
p6<-FeaturePlot(HCPUAB,features ="PTPRC",cols = c("#D3D3D3","#CD534C")) +
  labs(title = expression(Microglia~(italic("PTPRC"))),x = "UMAP1", y = "UMAP2")+
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 7))

pdf("output/summary/HCPUAB_celltype_FeaturePlot.pdf", height = 6, width =8 )
p<- (p1|p2|p3)/(p4|p5|p6)
dev.off()

library(patchwork)

figure <- (p1 | p2 | p3) / (p4 | p5 | p6)

ggsave(
  "output/summary/HCPUAB_celltype_FeaturePlot.pdf",
  plot = figure,
  width = 8,
  height = 6,
  device = cairo_pdf,
  family = "Arial"
)

DotPlot(HCPUAB,assay = 'RNA', features =c("AQP4","PLP1","PDGFRA","PTPRC","APOLD1"), dot.scale = 5) +
  RotatedAxis()
Maker.cellmarker2 <-c("AQP4","GFAP",#Astro
                      "PTPRC","CX3CR1", #Micro
                      "MYT1L","CUX2", #Neuron
                      "PDGFRA","PCDH15", #OPC
                      "PLP1","MBP","MOG", #Oligodendrocytes
                      "DCN","FLT1" #Vas
)
DotPlot(HCPUAB,assay = 'RNA', features =Maker.cellmarker2, dot.scale = 5) +
  RotatedAxis()

TOPDEG_marker<-c("AQP4","GFAP",#Astro
                 "PTPRC","CX3CR1", #Micro
                 "MYT1L","CUX2", #Neuron
                 "GPR17","ADAM33","SEMA3D","CFAP95","ANLN","SLCO1A2", #Oligodendrocytes
                 "PDGFRA","MYT1","EGFR","GAS1", #OPC
                 "DCN","EBF1" #Vas
)
DotPlot(HCPUAB,assay = 'RNA', features =TOPDEG_marker, dot.scale = 5) +
  RotatedAxis()


DimPlot(HCPUAB,label = T)
HCPUAB$sampleID<- factor(HCPUAB$SampleID, 
                         levels = c("control_m1", "HCPUAB023", "control_y10", "HCPUAB014"))
DimPlot(HCPUAB,label = T,split.by = "orig.ident",ncol=2)


Top.DEG<- DEGs %>% group_by(group1) %>% top_n(n = 3, wt = avg_log2FC)
DotPlot(HCPUAB,assay = 'RNA', features =Maker.cellmarker2, dot.scale = 5) +
  RotatedAxis()


pdf("output/summary/cellmarker.dot.pdf", height = 4, width=8)

Maker.cellmarker2 <-c("AQP4","GFAP",#Astro
                      "PTPRC","CX3CR1", #Micro
                      "MYT1L","CUX2", #Neuron
                      "PLP1","MBP","MOG", #ODC
                      "PDGFRA","PCDH15" #OPC
)

avg_expression = DotPlot(HCPUAB, features = unique(Maker.cellmarker2))
avg_expression = avg_expression$data
avg_expression$features.plot <- factor(avg_expression$features.plot, levels = (rownames(avg_expression)))
avg_expression$id <- factor(avg_expression$id, levels = c("Astrocytes","Microglia","Neurons","ODC","OPC"))
ggplot(avg_expression,
       aes(x = features.plot, y = id, size = pct.exp, fill = avg.exp.scaled)) + 
  #geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  geom_point(shape = 21) +
  #scale_size_continuous(range = c(1, 4)) +
  #scale_fill_gradientn(colours = pals::coolwarm(100)) +
  scale_fill_gradient2(low = "#053061",
                       mid = "#eeeeee",
                       high = "#67001F") +
  ylab("") +
  xlab("") +
  theme_minimal() +
  #coord_fixed(ratio=3) +
  labs(size="Percent", fill="Avg. Exp") +
  theme(axis.text = element_text(color="black"), text = element_text(size=10),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        legend.key.size = unit(0.3, 'cm'), legend.position="bottom")+
  coord_flip() #翻转

dev.off()

pdf("output/summary/cellmarker.violin.pdf", height = 6, width=6)
FeatureStatPlot(HCPUAB,slot = "counts", stat.by = Maker.cellmarker2, group.by = "celltype", bg.by = "celltype",  plot_type ="violin", stack = TRUE,palcolor =my_cols,bg_palcolor = my_cols)
dev.off()

save(HCPUAB, file="output/RData/3_integrate_HCPUAB.RData")
saveRDS(HCPUAB, file="output/RData/3_integrate_HCPUAB.rds")




setwd("~/Library/CloudStorage/OneDrive-UAB-TheUniversityofAlabamaatBirmingham/1_ChongLab/Project4_Andy_hydrocephalus/DuyPhan_genecheck/")

source(file = "scripts/Plot_colorPaletters.R")
source(file = "scripts/load_libraries.R")

load("output/RData/1_merged_unfiltered_seurat.RData")

metadata <- merged_seurat@meta.data


pdf(file = "output/summary/1_unfiltered_QC_1.pdf", width = 8, height = 6)
# Visualize the number of cell counts per sample
# Rename columns

metadata %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 400) +
  facet_wrap(~orig.ident,nrow = 5)

metadata %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)


# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 20)

VlnPlot(merged_seurat,features = 'nCount_RNA',pt.size = 0,group.by = 'orig.ident')+
  ylim(NA,30000)+geom_hline(yintercept = 200) +
  geom_boxplot(width=0.1, fill="white",pt.size = 0,outliers = F)
VlnPlot(merged_seurat,features = "nFeature_RNA",pt.size = 0,group.by = 'orig.ident')+
  ylim(NA,10000)+geom_hline(yintercept = 200) +
  geom_boxplot(width=0.1, fill="white",pt.size = 0,outliers = F)
VlnPlot(merged_seurat,features = "percent.mt",pt.size = 0,group.by = 'orig.ident')+
  ylim(NA,50)+geom_hline(yintercept = 20) +
  geom_boxplot(width=0.1, fill="white",pt.size = 0,outliers = F)

plot1 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")+NoLegend()
plot2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

dev.off()


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
pdf(file = "output/summary/1_unfiltered_QC_2.pdf", width = 14, height = 4)
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0, ncol = 3)
dev.off()

#### remove doublet
library(scDblFinder) # Require cleanup of low-quality cells in advance
#### split the dataset into a list of two seurat objects (stim and CTRL)
#### https://satijalab.org/seurat/archive/v3.1/immune_alignment.html


#CALCULATE DOUBLET SCORE
merged_seurat[["RNA"]] <- as(object = merged_seurat[["RNA"]], Class = "Assay")
doublet.score<- scDblFinder::scDblFinder(as.SingleCellExperiment(merged_seurat,assay = "RNA"), samples="orig.ident", BPPARAM=
                                           BiocParallel::MulticoreParam(3),returnType="table")
merged_seurat <- AddMetaData(object = merged_seurat, metadata = doublet.score[colnames(merged_seurat),"score"], col.name = "doublet.score")
merged_seurat <- AddMetaData(object = merged_seurat, metadata = doublet.score[colnames(merged_seurat),"class"], col.name = "doublet.class")
table(merged_seurat$doublet.class)
#Filter 3: remove doublet
table(merged_seurat$orig.ident)
filtered_seurat <- subset(merged_seurat, subset = doublet.class == "singlet")
filtered_seurat <- subset(filtered_seurat, subset = nFeature_RNA > 200 &nFeature_RNA< 7500& percent.mt < 5)

cell.num<- table(filtered_seurat$orig.ident)
write.csv(cell.num, file = "output/summary/filtered_cell_number.csv")
# Visualize the number of cell counts per sample
# Rename columns
metadata <- filtered_seurat@meta.data
pdf(file = "output/summary/2_filtered_QC_1.pdf", width = 8, height = 6)
# Visualize the number of cell counts per sample
# Rename columns

metadata %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 400) +
  facet_wrap(~orig.ident,nrow = 5)

metadata %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)


# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 20)

VlnPlot(merged_seurat,features = 'nCount_RNA',pt.size = 0,group.by = 'orig.ident')+
  ylim(NA,30000)+geom_hline(yintercept = 200) +
  geom_boxplot(width=0.1, fill="white",pt.size = 0,outliers = F)
VlnPlot(merged_seurat,features = "nFeature_RNA",pt.size = 0,group.by = 'orig.ident')+
  ylim(NA,10000)+geom_hline(yintercept = 200) +
  geom_boxplot(width=0.1, fill="white",pt.size = 0,outliers = F)
VlnPlot(merged_seurat,features = "percent.mt",pt.size = 0,group.by = 'orig.ident')+
  ylim(NA,50)+geom_hline(yintercept = 20) +
  geom_boxplot(width=0.1, fill="white",pt.size = 0,outliers = F)

plot1 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")+NoLegend()
plot2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

dev.off()


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
pdf(file = "output/summary/2_filtered_QC_2.pdf", width = 14, height = 4)
VlnPlot(filtered_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0, ncol = 3)
dev.off()
# Create .RData object to load at any time
save(filtered_seurat, file="output/RData/2_merged_filtered_seurat.RData")

#### normalize and identify variable features for each dataset independently
filtered_seurat@assays$RNA <-split(filtered_seurat@assays$RNA,f=filtered_seurat$orig.ident)
filtered_seurat <- NormalizeData(filtered_seurat)
filtered_seurat <- FindVariableFeatures(filtered_seurat)
filtered_seurat <- ScaleData(filtered_seurat)
filtered_seurat <- RunPCA(filtered_seurat)
#filtered_seurat <- FindNeighbors(filtered_seurat, dims = 1:30, reduction = "pca")
#filtered_seurat <- FindClusters(filtered_seurat, resolution = 0.4, cluster.name = "unintegrated_clusters")
#filtered_seurat <- RunUMAP(filtered_seurat, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
#DimPlot(filtered_seurat, reduction = "umap.unintegrated", group.by = c("orig.ident", "unintegrated_clusters"))

HCPUAB <- IntegrateLayers(
  object = filtered_seurat, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
HCPUAB  <- JoinLayers(HCPUAB)



#check cell number
cellnumber_sample_Etiology<-table(HCPUAB@meta.data$orig.ident)
write.csv(cellnumber_sample_Etiology,file="output/summary/cellnumber_sample_Etiology.csv")


#set dims 
ElbowPlot(HCPUAB,ndims = 50) #select 50
p1 <-DimPlot(HCPUAB, reduction = "pca", group.by = "orig.ident")
p2 <-DimPlot(HCPUAB, reduction = "harmony", group.by = "orig.ident")
p1|p2

HCPUAB <- FindNeighbors(HCPUAB, reduction = "harmony", dims = 1:50)
HCPUAB<- FindClusters(HCPUAB, resolution = c(0.1,0.2,0.3,0.4,0.6,1))


library(clustree)
pdf("output/summary/cluster_tree.pdf", width = 8, height = 8)
clustree(HCPUAB@meta.data, prefix = "RNA_snn_res.")
dev.off()

HCPUAB <- FindClusters(HCPUAB, resolution = 0.3)
table(HCPUAB$seurat_clusters)
HCPUAB <- RunUMAP(HCPUAB, dims = 1:30, reduction = "harmony")
DimPlot(HCPUAB, reduction = "umap", group.by = c("seurat_clusters"))

save(HCPUAB, file="output/RData/3_integrate_HCPUAB.RData")





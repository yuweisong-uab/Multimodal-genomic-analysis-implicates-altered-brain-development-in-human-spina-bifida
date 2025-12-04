library(cellchat)
setwd("~/Library/CloudStorage/OneDrive-UAB-TheUniversityofAlabamaatBirmingham/1_ChongLab/Project4_Andy_hydrocephalus/DuyPhan_genecheck/")
## ==== load 4 per-sample CellChat objects ====
cellchat.control_m1  <- readRDS("output/summary/cellchat/4_cellchat_control_m1.rds")
cellchat.control_y10 <- readRDS("output/summary/cellchat/4_cellchat_control_y10.rds")
cellchat.Patient_m1  <- readRDS("output/summary/cellchat/4_cellchat_HCPUAB023.rds")
cellchat.Patient_y10 <- readRDS("output/summary/cellchat/4_cellchat_HCPUAB014.rds")

## Optional: compute centrality for each if not already done
cellchat.control_m1  <- netAnalysis_computeCentrality(cellchat.control_m1,  slot.name = "netP")
cellchat.control_y10 <- netAnalysis_computeCentrality(cellchat.control_y10, slot.name = "netP")
cellchat.Patient_m1  <- netAnalysis_computeCentrality(cellchat.Patient_m1,  slot.name = "netP")
cellchat.Patient_y10 <- netAnalysis_computeCentrality(cellchat.Patient_y10, slot.name = "netP")

## ==== put them into a list ====
## 2) put into a list in the order you want:
object.list <- list(
  "1M control"  = cellchat.control_m1,
  "1M patient"  = cellchat.Patient_m1,
  "10Y control" = cellchat.control_y10,
  "10Y patient" = cellchat.Patient_y10
)
## ==== merge them ====
ptm <- Sys.time()
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
execution.time <- Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

## (optional) save merged
save(object.list, file = "output/summary/cellchat/cellchat_object_4samples.RData")
save(cellchat,   file = "output/summary/cellchat/cellchat_merged_4samples.RData")
## cell type colors
my_cols <- c(
  'Astrocytes'          = '#8F7700',
  'OPC'                 = '#3B3B3B',
  'Oligodendrocytes'    = '#EFC001',
  'Microglia'           = '#CD534C',
  'Excitatory neurons'  = '#59753E',
  'Inhibitory neurons'  = '#394E2F'
)


## 1M pair
object.list.1M  <- object.list[c("1M control", "1M patient")]
cellchat.1M <- mergeCellChat(object.list.1M, add.names = names(object.list.1M))

## 10Y pair
object.list.10Y <- object.list[c("10Y control", "10Y patient")]
cellchat.10Y <- mergeCellChat(object.list.10Y, add.names = names(object.list.10Y))

## barplots of interactions
gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = 1:4)
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = 1:4, measure = "weight")
gg1
ggsave("output/summary/cellchat/compareInteractions.pdf", width = 8, height = 5)
gg2
ggsave("output/summary/cellchat/InteractionsStrength.pdf", width = 8, height = 5)

## circle diff plot
pdf("output/summary/cellchat/1M_diffInteractions.pdf", width = 6, height = 5)
par(mfrow = c(1,1), xpd = TRUE)
netVisual_diffInteraction(cellchat.1M, weight.scale = TRUE, title.name=" Differential number of interactions \n 1M patient versus 1M control")
dev.off()
pdf("output/summary/cellchat/10Y_diffInteractions.pdf", width = 6, height = 5)
par(mfrow = c(1,1), xpd = TRUE)
netVisual_diffInteraction(cellchat.10Y, weight.scale = TRUE,title.name = " Differential number of interactions \n 10Y patient versus 10Y control")
dev.off()


## heatmap
h1_1M <- netVisual_heatmap(cellchat.1M)
h2_1M <- netVisual_heatmap(cellchat.1M, measure = "weight")
h1_1M + h2_1M

h1_10Y <- netVisual_heatmap(cellchat.10Y)
h2_10Y <- netVisual_heatmap(cellchat.10Y, measure = "weight")
h1_10Y + h2_10Y
ggsave("output/summary/cellchat/1M_compareInteractions.pdf", width = 8, height = 5)


# all sources/targets, dataset 1 vs 2
netVisual_bubble(
  cellchat.1M,
  sources.use = 1,
  targets.use = 1:6,
  comparison  = c(1, 2),
  angle.x     = 45
)
















pdf("output/summary/Myelomeningocele_CTRL_CirclePlot_4samples.pdf", width = 10, height = 10)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(2,2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$count,
    weight.scale     = TRUE,
    label.edge       = FALSE,
    edge.label.cex   = 0.6,
    edge.weight.max  = weight.max[2],
    edge.width.max   = 12,
    title.name       = names(object.list)[i],
    color.use        = my_cols
  )
}
dev.off()


cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

rankSimilarity(cellchat, type = "functional")











cellchat.control_m1<- readRDS("output/summary/cellchat/4_cellchat_control_m1.rds")
cellchat.control_y10<- readRDS("output/summary/cellchat/4_cellchat_control_y10.rds")
cellchat.Patient_m1<- readRDS("output/summary/cellchat/4_cellchat_HCPUAB023.rds")
cellchat.Patient_y10<- readRDS("output/summary/cellchat/4_cellchat_HCPUAB014.rds")


cellchat.Control   <- netAnalysis_computeCentrality(cellchat.Control, slot.name = "netP")



object.list <- list(Control = cellchat.Control,Patient = cellchat.Patient)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
# Users can now export the merged CellChat object and the list of the two separate objects for later use
#save(object.list, file = "cellchat_object.Myelomeningocele_CTRL.RData")
#save(cellchat, file = "cellchat_merged_Myelomeningocele_CTRL.RData")
my_cols <- c('Astrocytes'='#8F7700','OPC'='#3B3B3B','Oligodendrocytes'='#EFC001','Microglia'='#CD534C','Excitatory neurons'='#59753E','Inhibitory neurons'='#394E2F')
#my_cols <- c("Astrocytes"="#CD534C", "Microglia"="#8F7700","Oligodendrocytes"="#EFC001","OPC"="#93ADC0",'Excitatory neurons'='#59753E','Inhibitory neurons'='#2F4227')
load("cellchat_merged_Myelomeningocele_CTRL.RData")
ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#Circle plot showing differential number of interactions or interaction strength among different cell populations across two datasets
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

#Circle plot showing the number of interactions or interaction strength among different cell populations across multiple datasets
pdf("output/summary/Myelomeningocele_CTRL_CirclePlot.pdf", width = 10, height = 5)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F,edge.label.cex = 0.6, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0(names(object.list)[i]),color.use=my_cols)
}
dev.off()


ptm = Sys.time()
#Identify signaling groups based on their functional similarity

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

#netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)
#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
#Compute and visualize the pathway distance in the learned joint manifold
rankSimilarity(cellchat, type = "functional")
#> Compute the distance of signaling networks between datasets 1 2

#Identify altered signaling with distinct interaction strength
pdf("output/summary/cellchat_top_signaling_pathways.pdf", width = 6, height = 10)
rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE, thresh = 0.05)
#rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
dev.off()

#Compare outgoing (or incoming) signaling patterns associated with each cell population
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 18,color.use=my_cols)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 18,color.use=my_cols)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

pdf("output/summary/cellchat_signalingRole_heatmap.pdf", width = 8, height = 10)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 18, color.heatmap = "OrRd",angle.x = 45)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 18, color.heatmap = "OrRd",angle.x = 45)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()


#Identify dysfunctional signaling by comparing the communication probabities
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:6),  comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

#Identify dysfunctional signaling by using differential expression analysis
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Patient"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "Patient",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "Control",ligand.logFC = -0.05, receptor.logFC = NULL)
net.up.top50 <- net.up[order(net.up$ligand.logFC, decreasing = TRUE), ][1:50, ]
net.down.top50 <- net.down[order(net.down$ligand.logFC, decreasing = FALSE), ][1:50, ]



gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
#df <- findEnrichedSignaling(object.list[[2]], features = c("CCL19", "CXCL12"), idents = c("Inflam. FIB", "COL11A1+ FIB"), pattern ="outgoing")


# Chord diagram
pdf("output/summary/cellchat_chord_differential_signaling_up.pdf", width = 6, height = 6)
netVisual_chord_gene(object.list[[2]],  slot.name = 'net', net = net.up.top50, lab.cex = 0.6, small.gap = 0.5,reduce=0.005,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]),color.use=my_cols,show.legend=F)
#> You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway
dev.off()

pdf("output/summary/cellchat_chord_differential_signaling_down.pdf", width = 6, height = 6)
netVisual_chord_gene(object.list[[1]],  slot.name = 'net', net = net.down.top50, lab.cex = 0.6, small.gap = 0.5,reduce=0.005, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]),color.use=my_cols,show.legend=F)
dev.off()


netAnalysis_river(object.list[[2]], pattern = "outgoing")

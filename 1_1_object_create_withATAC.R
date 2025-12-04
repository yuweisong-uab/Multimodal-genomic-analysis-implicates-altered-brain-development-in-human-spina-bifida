setwd("~/Library/CloudStorage/OneDrive-UAB-TheUniversityofAlabamaatBirmingham/1_ChongLab/Project4_Andy_hydrocephalus/DuyPhan_genecheck/")

source(file = "scripts/Plot_colorPaletters.R")
source(file = "scripts/setupseurat_function.R")
source(file = "scripts/load_libraries.R")
options(Seurat.object.assay.version = "v3")
#generate all seurat objects
lapply(c("HCPUAB014"),setupseurat)
lapply(c("HCPUAB023","control_m1","control_y10"),setupseurat_GEX)
setwd("~/Library/CloudStorage/OneDrive-UAB-TheUniversityofAlabamaatBirmingham/1_ChongLab/Project4_Andy_hydrocephalus/scRNA_analysis/")
# Create a Seurat object for each sample
for (file in c("HCPUAB014","HCPUAB023","control_m1","control_y10")){
  seurat_obj <-readRDS(paste0("cell_ranger_output/", file,"/",file,".SeuratObject.rds"))
  assign(file, seurat_obj)
}

#merge data
merged_seurat <- merge(x = HCPUAB014, 
                       y = c(HCPUAB023,control_m1,control_y10),
                       add.cell.id = c("HCPUAB014","HCPUAB023","control_m1","control_y10"))

setwd("~/Library/CloudStorage/OneDrive-UAB-TheUniversityofAlabamaatBirmingham/1_ChongLab/Project4_Andy_hydrocephalus/DuyPhan_genecheck/")
# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)
cell.num<- table(merged_seurat$orig.ident)
write.csv(cell.num, file = "output/initial_cell_number.csv")
# Explore merged metadata
View(merged_seurat@meta.data)
# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
# Compute percent mito ratio
merged_seurat$percent.mt <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$percent.ribo<- PercentageFeatureSet(object = merged_seurat, pattern = "^RP[SL]")
# Create metadata dataframe
metadata <- merged_seurat@meta.data


# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
# Create sample column
metadata$Etiology <- NA
metadata$Etiology[which(str_detect(metadata$cells, "^HCPUAB014_"))] <- "Myelomeningocele"
metadata$Etiology[which(str_detect(metadata$cells, "^HCPUAB023_"))] <- "Myelomeningocele"
metadata$Etiology[which(str_detect(metadata$cells, "^control_m1_"))] <- "CTRL_m1"
metadata$Etiology[which(str_detect(metadata$cells, "^control_y10_"))] <- "CTRL_y10"

metadata$Race <- NA
metadata$Race[which(str_detect(metadata$cells, "^HCPUAB014_"))] <- "Afraican American"
metadata$Race[which(str_detect(metadata$cells, "^HCPUAB023_"))] <- "Caucasian"
metadata$Race[which(str_detect(metadata$cells, "^control_m1_"))] <- "Caucasian"
metadata$Race[which(str_detect(metadata$cells, "^control_y10_"))] <- "Afraican American"

metadata$Gender <- NA
metadata$Gender[which(str_detect(metadata$cells, "^HCPUAB014_"))] <- "Male"
metadata$Gender[which(str_detect(metadata$cells, "^HCPUAB023_"))] <- "Female"
metadata$Gender[which(str_detect(metadata$cells, "^control_m1_"))] <- "Female"
metadata$Gender[which(str_detect(metadata$cells, "^control_y10_"))] <- "Male"

metadata$Ages <- NA
metadata$Ages[which(str_detect(metadata$cells, "^HCPUAB014_"))] <- "y10"
metadata$Ages[which(str_detect(metadata$cells, "^HCPUAB023_"))] <- "m1"
metadata$Ages[which(str_detect(metadata$cells, "^control_m1_"))] <- "m1"
metadata$Ages[which(str_detect(metadata$cells, "^control_y10_"))] <- "y10"


metadata$SampleID <- NA
metadata$SampleID[which(str_detect(metadata$cells, "^HCPUAB014_"))] <- "HCPUAB014"
metadata$SampleID[which(str_detect(metadata$cells, "^HCPUAB023_"))] <- "HCPUAB023"
metadata$SampleID[which(str_detect(metadata$cells, "^control_m1_"))] <- "control_m1"
metadata$SampleID[which(str_detect(metadata$cells, "^control_y10_"))] <- "control_y10"

metadata$Condition <- NA
metadata$Condition[which(str_detect(metadata$cells, "^HCPUAB014_"))] <- "Patient"
metadata$Condition[which(str_detect(metadata$cells, "^HCPUAB023_"))] <- "Patient"
metadata$Condition[which(str_detect(metadata$cells, "^control_m1_"))] <- "Control"
metadata$Condition[which(str_detect(metadata$cells, "^control_y10_"))] <- "Control"
metadata$SampleID <- NA
metadata$SampleID[which(str_detect(metadata$cells, "^HCPUAB014_"))] <- "10Y patient"
metadata$SampleID[which(str_detect(metadata$cells, "^HCPUAB023_"))] <- "1M patient"
metadata$SampleID[which(str_detect(metadata$cells, "^control_m1_"))] <- "10Y control"
metadata$SampleID[which(str_detect(metadata$cells, "^control_y10_"))] <- "1M control"

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

head(merged_seurat@meta.data)

# Create .RData object to load at any time
save(merged_seurat, file="output/RData/1_merged_unfiltered_seurat.RData")


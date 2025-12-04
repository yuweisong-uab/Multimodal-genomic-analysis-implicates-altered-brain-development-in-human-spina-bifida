## Load Seurat and dplyr
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
## Setup parallelization (https://cran.r-project.org/web/packages/future/index.html)
library(future)
plan("multicore", workers=4)
library(CellChat)

#renv::activate(project = "/Users/yuweisong/miniconda3/envs/SCP_env")
plan("multisession", workers =8) 
options(future.globals.maxSize = 50000 * 1024^2) # set 50G RAM
setwd("~/Library/CloudStorage/OneDrive-UAB-TheUniversityofAlabamaatBirmingham/1_ChongLab/Project4_Andy_hydrocephalus/DuyPhan_genecheck/")
load("output/RData/3_integrate_HCPUAB.RData")

Idents(HCPUAB)<-HCPUAB$celltype
table(HCPUAB$SampleID)
table(HCPUAB$celltype)

HCPUAB[["RNA3"]] <- as(object = HCPUAB[["RNA"]], Class = "Assay")
DefaultAssay(HCPUAB) <- "RNA3"
HCPUAB[["RNA"]] <- NULL
HCPUAB <- RenameAssays(object = HCPUAB, RNA3 = 'RNA')



control_m1_rds <- subset(HCPUAB, subset = SampleID %in% c('control_m1'))
table(control_m1_rds$SampleID)

ptm = Sys.time()
data.input <- control_m1_rds[["RNA"]]@data # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
labels <- Idents(control_m1_rds)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
#Create a CellChat object
cellchat <- createCellChat(object = control_m1_rds, group.by = "celltype", assay = "RNA")
#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use
#Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
#Part II: Inference of cell-cell communication network
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
# par(mfrow = c(3,4), xpd=TRUE)
# for (i in 1:nrow(mat)) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
# }

saveRDS(cellchat, file = "output/RData/4_cellchat_control_m1.rds")




run_cellchat_for_sample <- function(
    seurat_obj,
    sample_id,
    group_by   = "celltype",
    assay      = "RNA",
    species    = c("human", "mouse"),
    workers    = 4,
    min_cells  = 10,
    out_dir    = "output/RData",
    do_plots   = TRUE
) {
  species <- match.arg(species)
  
  # ---- checks ----
  if (!inherits(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.")
  }
  if (!"SampleID" %in% colnames(seurat_obj@meta.data)) {
    stop("`SampleID` not found in seurat_obj@meta.data.")
  }
  if (!sample_id %in% seurat_obj$SampleID) {
    stop(sprintf("`sample_id` ('%s') not found in seurat_obj$SampleID. Available: %s",
                 sample_id, paste(sort(unique(seurat_obj$SampleID)), collapse = ", ")))
  }
  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("`group_by` ('%s') not found in seurat_obj@meta.data.", group_by))
  }
  if (!assay %in% names(seurat_obj@assays)) {
    stop(sprintf("Assay '%s' not found in seurat_obj@assays.", assay))
  }
  
  # ---- subset to the requested sample ----
  message("Subsetting to SampleID == ", sample_id)
  obj <- subset(seurat_obj, subset = SampleID %in% sample_id)
  message("Cell counts by SampleID after subset:")
  print(table(obj$SampleID))
  
  # ---- choose database ----
  suppressPackageStartupMessages(require(CellChat))
  CellChatDB <- if (species == "human") CellChatDB.human else CellChatDB.mouse
  CellChatDB.use <- CellChatDB
  
  # ---- create CellChat object ----
  message("Creating CellChat object …")
  cellchat <- createCellChat(object = obj, group.by = group_by, assay = assay)
  
  # set DB
  cellchat@DB <- CellChatDB.use
  showDatabaseCategory(CellChatDB.use)
  
  # ---- timing helpers ----
  t0 <- Sys.time()
  tick <- function(label) {
    dt <- as.numeric(Sys.time() - t0, units = "secs")
    message(sprintf("[%s] elapsed: %.2f sec", label, dt))
  }
  
  # ---- parallel plan (restore on exit) ----
  old_plan <- future::plan()
  on.exit({
    try(future::plan(old_plan), silent = TRUE)
  }, add = TRUE)
  future::plan("multisession", workers = workers)
  
  # ---- Part I: preprocessing ----
  message("Preprocessing …")
  cellchat <- subsetData(cellchat)                      # subset signaling genes
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  tick("Preprocessing done")
  
  # ---- Part II: inference ----
  message("Inferring communication probabilities and pathways …")
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  tick("Inference done")
  
  # ---- plots (optional) ----
  if (isTRUE(do_plots)) {
    message("Plotting circle networks …")
    groupSize <- as.numeric(table(cellchat@idents))
    op <- par(mfrow = c(1,2), xpd = TRUE)
    on.exit(par(op), add = TRUE)
    netVisual_circle(cellchat@net$count,
                     vertex.weight = groupSize,
                     weight.scale  = TRUE,
                     label.edge    = FALSE,
                     title.name    = "Number of interactions")
    netVisual_circle(cellchat@net$weight,
                     vertex.weight = groupSize,
                     weight.scale  = TRUE,
                     label.edge    = FALSE,
                     title.name    = "Interaction weights/strength")
  }
  
  # ---- save RDS ----
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_file <- file.path(out_dir, sprintf("4_cellchat_%s.rds", sample_id))
  saveRDS(cellchat, file = out_file)
  message("Saved: ", out_file)
  
  invisible(cellchat)
}

# Vector of sample IDs you want to run
#"control_m1",
samples_to_run <- c( "control_y10", "HCPUAB014", "HCPUAB023")

# Loop over them
cellchat_list <- list()
for (sid in samples_to_run) {
  message("\n==== Running CellChat for ", sid, " ====")
  cellchat_list[[sid]] <- run_cellchat_for_sample(
    seurat_obj = HCPUAB,
    sample_id  = sid,
    group_by   = "celltype",
    assay      = "RNA",
    species    = "human",     # change to "mouse" if needed
    workers    = 4,
    min_cells  = 10,
    out_dir    = "output/RData",
    do_plots   = TRUE
  )
}


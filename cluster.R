#'10x data cluster (intergrate different datasets)
#'Data: 2019/11/18
#'Author: Han Zhifa

#*####**************************************####*#
#             packages and workDir               #
######======================================######
#.libPaths()
#install.packages("BiocManager")
#BiocManager::install("scater")
#install.packages("mvoutlier")
#install.packages("cellranger")
#require(cellranger)
#require(Matrix)
#("pheatmap")
#install.packages("clustree")
try(library(Seurat, lib.loc = "/home/hanzf/.conda/envs/scs10x/lib/R/library"), 
    silent = T)
require(Seurat, lib.loc = "/opt/R/4.1.0/lib/R/library")

require(cowplot)

#getwd()
setwd("/ssd/hanzf/research/scs10x")

#####

#*####**************************************####*#
#                    get data                    #
######======================================######
## get results outputing path
results_path <- commandArgs(T)[1]

## read combined
combined <- readRDS(commandArgs(T)[2])

## read PCs
PCs <- read.csv(commandArgs(T)[3])$PCs

## get resolution of cluster
resolution <- read.csv(commandArgs(T)[4])$resolution

#####

#*####**************************************####*#
#     DimPlot for different reduction method     #
######======================================######
if("integrated" %in% names(combined@assays) ){
  DefaultAssay(combined) <- "integrated"
}

combined <- FindNeighbors(combined, reduction = "pca", dims = PCs)
combined <- FindClusters(combined, resolution = resolution)
## plot PCA, tSNE and UMAP
# run PCA
combined <- RunPCA(combined, npcs = PCs[length(PCs)])
# plot PCA
p1 <- DimPlot(object = combined, reduction = "pca", group.by = "experiment")
p2 <- DimPlot(combined, reduction = "pca", label = TRUE)
png(paste0(results_path, "PCA_cluster.png"), width = 1300)
print(plot_grid(p1, p2))
dev.off()
# run tSNE
combined <- RunTSNE(combined, dims.use = PCs, check_duplicates = FALSE)
# plot tSNE
p1 <- DimPlot(object = combined, reduction = "tsne", group.by = "experiment")
p2 <- DimPlot(combined, reduction = "tsne", label = TRUE)
png(paste0(results_path, "tSNE_cluster.png"), width = 1300)
print(plot_grid(p1, p2))
dev.off()
# run UMAP
combined <- RunUMAP(combined, reduction = "pca", dims = PCs)
# Visualization UMAP
p1 <- DimPlot(combined, reduction = "umap", group.by = "experiment")
p2 <- DimPlot(combined, reduction = "umap", label = TRUE)
png(paste0(results_path, "UMAP_cluster.png"), width = 1300)
print(plot_grid(p1, p2))
dev.off()
# plot UMAP split by experiment
if(length(unique(combined$experiment)) > 3){
  png(paste0(results_path, "UMAP_splitByExperiment.png"), 
      height = length(unique(combined$seurat_clusters)) /3 * 90, 
      width = 1500)
  print(DimPlot(combined, reduction = "umap", split.by = "experiment", ncol = 3))
  dev.off()
}else{
  png(paste0(results_path, "UMAP_splitByExperiment.png"))
  print(DimPlot(combined, reduction = "umap", split.by = "experiment"))
  dev.off()
}
# plot UMAP split by cluster
png(paste0(results_path, "UMAP_splitByCluster.png"), 
    height = length(unique(combined$seurat_clusters))/3 * 170, width = 900)
print(DimPlot(combined, reduction = "umap", split.by = "seurat_clusters", ncol = 3))
dev.off()

#####

#*####**************************************####*#
#       plot doublet rank for each cluster       #
######======================================######
pdf(paste0(results_path, "doublet_rank.pdf"), 
    width = 3,#length(unique(combined$seurat_clusters)) /3 * 1.7, 
    height = 7)
VlnPlot(combined, features = c("doublet_scores"), 
        split.by = "seurat_clusters")
dev.off()

#####

## save rds file
#file.remove(commandArgs(T)[2])
saveRDS(combined, file = paste0(results_path, "combined_cluster.rds"))

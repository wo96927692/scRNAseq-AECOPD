#'10x data cluster (intergrate different datasets)
#'Data: 2021/9/26
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

#getwd()
setwd("/ssd/hanzf/research/scs10x") 

#####

#*####**************************************####*#
#                    get data                    #
######======================================######
## get results outputing path
results_path <- commandArgs(T)[1]
LT_path <- paste0(results_path, "LabelTransfer/")
dir.create(LT_path)

## query seurat data
combined <- readRDS(file = commandArgs(T)[2])

## reference rds file name
cat(paste0("Whether the reference seurat object stored in ", LT_path, "?(T or F):"))
if(as.logical(readLines("stdin", n=1))){
  cat("input the rds file name of the seurat class which will be label transfer reference: \n")
  file_name <- readLines("stdin", n = 1)
  reference_seurat <- readRDS(file = paste0(LT_path, file_name))
}else{
  cat("input the rds file stored path and file name of the seurat class which will be label transfer reference: \n")
  path_name <- readLines("stdin", n = 1)
  reference_seurat <- readRDS(file = path_name)
}

## PCs
PCs <- read.csv(commandArgs(T)[3], header = T, quote = "\"")$PCs

#####
if("integrated" %in% names(reference_seurat@assays) ){
  DefaultAssay(reference_seurat) <- "integrated"
}else{
  DefaultAssay(reference_seurat) <- "RNA"
}
if("integrated" %in% names(combined@assays) ){
  DefaultAssay(combined) <- "integrated"
}else{
  DefaultAssay(combined) <- "RNA"
}

## transfer
TF.anchors <- FindTransferAnchors(reference = reference_seurat, query = combined, 
                                        dims = PCs)
predictions <- TransferData(anchorset = TF.anchors, refdata = reference_seurat$celltypes, 
                            dims = PCs)
combined <- AddMetaData(combined, metadata = 
                          predictions[c("predicted.id", "prediction.score.max")])

## have a look
print(table(combined$predicted.id))

## PCA
# PCA
png(paste0(LT_path, "PCA.png"), 
    width = 480*2, height = 480*2)
print(DimPlot(combined, reduction = "pca", group.by = "predicted.id", label = T))
dev.off()
# tSNE
png(paste0(LT_path, "tSNE.png"), 
    width = 480*2, height = 480*2)
print(DimPlot(combined, reduction = "tsne", group.by = "predicted.id", label = T))
dev.off()
# UMAP
png(paste0(LT_path, "UMAP.png"), 
    width = 480*2, height = 480*2)
print(DimPlot(combined, reduction = "umap", group.by = "predicted.id", label = T))
dev.off()


#####
## savedata
saveRDS(combined, file = paste0(results_path, "combined_cluster.rds"))
write.csv(predictions, file = paste0(LT_path, "Label_predictions.csv"))

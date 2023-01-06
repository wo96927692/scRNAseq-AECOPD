#'10x data celltypes
#'Data: 2021/09/28
#'Author: Han Zhifa

#*####**************************************####*#
#             packages and workDir               #
######======================================######
#.libPaths()
require(Seurat)
require(ggplot2)
#library(dplyr)

#getwd()
#setwd("/ssd/hanzf/research/scs10x")

#####

#*####**************************************####*#
#                    get data                    #
######======================================######
##
results_path <- commandArgs(T)[2]
celltype_path <- paste0(results_path, "celltypes/")
dir.create(celltype_path)

## read combined
combined <- readRDS(commandArgs(T)[1])
if("integrated" %in% names(combined@assays) ){
  DefaultAssay(combined) <- "integrated"
}
Idents(combined) <- combined$seurat_clusters

#####

#*####**************************************####*#
#       get the celltypes for each cluster       #
######======================================######
## get the celltypes according clusters
cat(paste0("Whether the correct celltypes.csv file stored in ", celltype_path, " ?(T or F):"))
wthCellTypeFile <- as.logical(readLines("stdin", n=1))
if(wthCellTypeFile){
  celltypes_df <- read.csv(file = paste0(celltype_path, "celltypes.csv"), header = T, 
                           stringsAsFactors = F)
  clusters <- as.character(celltypes_df$clusters)
  celltypes <- celltypes_df$celltypes
  if(identical(clusters, levels(combined))){
    names(celltypes) <- clusters
  }else{
    cat("the levels of combined are:\n")
    cat(levels(combined))
    
    ## get the main celltypes for each cluster
    clusters <- levels(combined)
    celltypes <- character(length(clusters))
    names(celltypes) <- clusters
    for(i in clusters) {
      cat(paste0("cluster ", i, " is: "))
      celltypes[i] <- readLines("stdin", n = 1)
    }
  }
}else{
  ## get the main celltypes for each cluster
  clusters <- levels(combined)
  celltypes <- character(length(clusters))
  names(celltypes) <- clusters
  for(i in clusters) {
    cat(paste0("cluster ", i, " is: "))
    celltypes[i] <- readLines("stdin", n = 1)
  }
}
cat(celltypes)

## add the celltypes to combined
combined <- RenameIdents(combined, celltypes)
# add to metadata
combined<- AddMetaData(combined, data.frame(celltypes = combined@active.ident))
# deleta the DEL
combined <- combined[, !("DEL" == Idents(combined))]

#####

#*####**************************************####*#
#                  visualization                  #
######======================================######
## PCA
pdf(paste0(celltype_path, "PCA_celltypes.png"), width = 9)
DimPlot(combined, reduction = "pca", label = TRUE)
dev.off()

## tSNE
pdf(paste0(celltype_path, "tSNE_celltypes.png"), width = 9)
DimPlot(combined, reduction = "tsne", label = TRUE)
dev.off()

## UMAP
pdf(paste0(celltype_path, "UMAP_celltypes.png"), width = 9)
DimPlot(combined, reduction = "umap", label = TRUE)
dev.off()

## barplot for celltype percentage
table(combined$experiment)
png(paste0(celltype_path, "barplot.png"))
ggplot(combined@meta.data[, c("celltypes", "experiment")], 
       aes(x = factor(experiment, levels = unique(experiment)))) + 
  geom_bar(aes(fill = celltypes, 
               color = 1), width = .55, position = "fill") + 
  labs(x = "experiment")
dev.off()

## heatmap for celltype percentage

## save RDS
saveRDS(combined, file = paste0(results_path, "combined_celltypes.rds"))

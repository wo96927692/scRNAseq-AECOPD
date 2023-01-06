#'10x data cluster (intergrate different datasets)
#'Data: 2021/09/26
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
#install.packages("pheatmap")
try(library(Seurat, lib.loc = "/home/hanzf/.conda/envs/scs10x/lib/R/library"), 
    silent = T)
require(Seurat, lib.loc = "/opt/R/4.1.0/lib/R/library")
require(pheatmap)
#require(cowplot)
#require(future)
library(dplyr, lib.loc = "/opt/R/4.1.0/lib/R/library")

#getwd()
setwd("/ssd/hanzf/research/scs10x")

## set parallelization
#cat("set the workers used in finding markers (1 worker need 1Gb memory for 10Mb seurat class): ")
#workers <- as.numeric(readLines("stdin", n = 1))
#plan("multiprocess", workers = workers)
#options(future.globals.maxSize = workers * 1024^3 * 50)
# 
#plan()

#####

#*####**************************************####*#
#                    get data                    #
######======================================######
## 
results_path <- commandArgs(T)[2]

## read main celltypes
# function: read dictionary args like: "celltype1:marker1,melltype2:marker2"
readDic <- function(arg) {
  arg <- as.character(unlist(strsplit(arg, ",")))
  index <- 1:length(arg)
  List_second_split <- tapply(arg, index, function(x) {
    strsplit(x, ":")
  })
  names <- as.character()
  result <- as.character()
  for( i in index ) {
    names[i] <- List_second_split[[i]][1]
    result[i] <- List_second_split[[i]][2]
  }
  names(result) <- names
  return(result)
}
# read main cell
markers_celltypes_main <- readDic(commandArgs(T)[1])

## get resolution of cluster
resolution <- read.csv(commandArgs(T)[4])$resolution

## read PCs
PCs <- read.csv(commandArgs(T)[5])$PCs

## read combined
combined <- readRDS(commandArgs(T)[3])

## create marker path
marker_path <- paste0(results_path, "markers/")
dir.create(marker_path)

#####

#*####**************************************####*#
#             Plot for main markers              #
######======================================######
## identify markers for each cluster
DefaultAssay(combined) <- "RNA"

# FeaturePlot PCA
pdf(paste0(marker_path, "FeaturePlot_PCA_main_markers.pdf"))
print(FeaturePlot(combined, features = markers_celltypes_main, 
                  reduction = "pca"))
dev.off()
# FeaturePlot tSNE
pdf(paste0(marker_path, "FeaturePlot_tSNE_main_markers.pdf"))
print(FeaturePlot(combined, features = markers_celltypes_main, 
                  reduction = "tsne"))
dev.off()
# FeaturePlot UMAP
pdf(paste0(marker_path, "FeaturePlot_UMAP_main_markers.pdf"))
print(FeaturePlot(combined, features = markers_celltypes_main, 
                  reduction = "umap"))
dev.off()

## RidgePlot
pdf(paste0(marker_path, "RidegPlot_mainMarkers.pdf"), 
    height = length(unique(combined$seurat_clusters)) * 0.15 * 
      length(markers_celltypes_main) + 0.7, 
    width = length(markers_celltypes_main) * 3.9)
print(RidgePlot(combined, features = markers_celltypes_main))
dev.off()

## DotPlot
pdf(paste0(marker_path, "DotPlot_mainMarkers.pdf"), 
    width = length(markers_celltypes_main)* 0.3 + 5, 
    height = length(unique(combined$seurat_clusters)) * 0.19 + 3)
print(DotPlot(combined, features = rev(markers_celltypes_main), 
              group.by = "seurat_clusters") + RotatedAxis())
dev.off()

## VlnPlot
jpeg(paste0(marker_path, "VlnPlot_mainMarkers.jpeg"), quality = 100, 
     width = length(unique(combined$seurat_clusters)) * 45 + 200, 
     height = ceiling(length(markers_celltypes_main) / 2) * 300)
print(VlnPlot(object = combined, 
              features = markers_celltypes_main, 
              ncol = 2, pt.size = 0))
dev.off()
# pdf
pdf(paste0(marker_path, "VlnPlot_mainMarkers.pdf"), 
    height = length(markers_celltypes_main) / 3 * 7, 
    width = length(unique(combined$seurat_clusters)) * 0.3 * 3)
print(VlnPlot(object = combined, features = markers_celltypes_main, pt.size = 0))
dev.off()

## average expression for each main marker
#average_expression_main_markers <- 
#  AverageExpression(combined, features = markers_celltypes_main)$RNA
#print(class(average_expression_main_markers))
#write.csv(average_expression_main_markers, 
#          file = paste0(marker_path, "average_expression_main_markers.csv"), 
#          row.names = T)

#####

#*####**************************************####*#
#        main markers and cluster markers        #
######======================================######
## find markers
if(F){
  conserved_markers <- list()
  dir.create(paste0(marker_path, "conserved_markers/"))
  for(i in c(0:length(levels(Idents(combined))))){
    conserved_markers[[i + 1]] <- FindConservedMarkers(combined, ident.1 = i, 
                                                       grouping.var = "orig.ident")
    names(conserved_markers)[i + 1] <- (paste0("markers", i))
    write.csv(conserved_markers[[i + 1]], 
              file = paste0(marker_path, "conserved_markers/cluster", i, ".csv"), 
              row.names = T)
  }
}else{
  ## find all markers
  all_markers <- FindAllMarkers(object = combined, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(all_markers, file = paste0(marker_path, "all_markers.csv"), row.names = T)
  ## find all positive markers
  all_pos_markers <- FindAllMarkers(object = combined, only.pos = T, 
                                    min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(all_pos_markers, 
            file = paste0(marker_path, "all_pos_markers.csv"), 
            row.names = T)
  
  ## scaling and regressing technology effects
  #if(length(unique(combined$batch)) > 1){
  #  combined <- ScaleData(combined, 
  #                        vars.to.regress = c("nCounts_RNA", "percent.mt", "batch"))
  #}else{
  #  combined <- ScaleData(combined, vars.to.regress = c("nCounts_RNA", "percent.mt"))
  #}
  
  ## plot heatmap
  #top10 <- all_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
  #pdf(paste0(marker_path, "Heatmap_top10_markers.pdf"), 
  #    width = ncol(combined) * 0.002 + 5)
  #print(DoHeatmap(object = combined, features = top10$gene, label = TRUE))
  #dev.off()
}

# get the cluster markers about main celltypes
cluster_main_markers <- merge(data.frame(gene = markers_celltypes_main), 
                              all_markers, by = "gene")
write.table(cluster_main_markers, file = paste0(marker_path, "cluster_main_markers"), 
            sep = "\t", row.names = F)

#####

## save data
#saveRDS(combined, file = paste0(results_path, "combined_markers.rds"))
#save.image(paste0(results_path, "find_main_markers.RData"))
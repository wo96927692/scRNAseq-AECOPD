#'10x data subClusters
#'Data: 2019/11/18
#'Author: Han Zhifa

#*####**************************************####*#
#             packages and workDir               #
######======================================######
#.libPaths()
try(library(Seurat, lib.loc = "/home/hanzf/.conda/envs/scs10x/lib/R/library"), 
    silent = T)
require(Seurat, lib.loc = "/opt/R/4.1.0/lib/R/library")
#require(ggplot2)
#library(dplyr)

#getwd()
setwd("/ssd/hanzf/research/scs10x")

#####

#*####**************************************####*#
#                    get data                    #
######======================================######
## read combined
combined <- readRDS(commandArgs(T)[1])

## get and create output path
sub_celltype <- commandArgs(T)[2]
results_path <- "results/"
celltype_subpath <- paste0(results_path, "cluster_", sub_celltype)
dir.create(celltype_subpath)

#####

#*####**************************************####*#
#              get the sub clusters              #
######======================================######
## get the sub clusters for sub_celltype
cat(paste0(sub_celltype, " includes the clusters (sperate by ',' (eg: 1,2,3) or use ':' operation (eg: 1:3)): \n"))
sub_clusters <- readLines("stdin", n = 1)
if(grepl(":", sub_clusters)){
  sub_clusters <- eval(parse(text = sub_clusters))
}else{
  sub_clusters <- as.numeric(unlist(strsplit(sub_clusters, ",")))
}
cat(paste0("The ", sub_celltype, " contain: ")); cat(sub_clusters); cat("\n")

## 
combined_sub <- combined[, combined$seurat_clusters %in% sub_clusters]

##
if("integrated" %in% names(combined_sub@assays) ){
  combined_sub@assays$integrated <- NULL
  DefaultAssay(combined_sub) <- "RNA"}


## save data
saveRDS(combined_sub, file = commandArgs(T)[3])

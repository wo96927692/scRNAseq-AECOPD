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
#require(future)
try(library(Seurat, lib.loc = "/home/hanzf/.conda/envs/scs10x/lib/R/library"), 
    silent = T)
require(Seurat, lib.loc = "/opt/R/4.1.0/lib/R/library")
require(clustree)

#getwd()
setwd("/ssd/hanzf/research/scs10x")

## set parallelization
#cat("set the workers used in clustree: ")
#workers <- as.numeric(readLines("stdin", n = 1))
#plan("multicore", workers = workers)
#options(future.globals.maxSize = workers * 1024^3)
# 
#plan()

#####

#*####**************************************####*#
#                    get data                    #
######======================================######
## read results output path
results_path <- commandArgs(T)[2]

## read PCs
PCs <- as.numeric(read.csv(commandArgs(T)[3])$PCs)

## input rds file name
file_path_name <- commandArgs(T)[1]
# read data
combined <- readRDS(file_path_name)

##
if("integrated" %in% names(combined@assays) ){
  DefaultAssay(combined) <- "integrated"
  snn_term <- "integrated"
}else{
  snn_term <- "RNA"
  DefaultAssay(combined) <- "RNA"
}
print(DefaultAssay(combined))
#####

#*####**************************************####*#
#                    clustree                    #
######======================================######
## Clustering
# KNN
combined <- FindNeighbors(combined, reduction = "pca", dims = PCs)
# clustering
res_seq <- c(seq(0, 0.15, 0.05), seq(0.2, 1, 0.1), seq(1.5, 2, 0.5))
for(i in res_seq){
  print(paste0("res_seq: ", i))
  combined <- FindClusters(combined, resolution = i)
  print(colnames(combined@meta.data) )
}
# plot the clustree
pdf(paste0(results_path, "clustree.pdf"), height = length(res_seq) * 0.7, 
    width = 12)
clustree(combined, prefix = paste0(snn_term, "_snn_res."))
dev.off()

#####

#*####**************************************####*#
#                get resolution                  #
######======================================######
## get resolution of cluster
cat("The resolution used in clustering is: \n")
chares <- readLines("stdin", n = 1)
resolution <- as.numeric(chares)
cat("The used resolutions is: "); cat(resolution); cat("\n")

#####

## save RDS
write.csv(data.frame(resolution = resolution), file = paste0(results_path, "resolution.csv"), 
          row.names = F)
#saveRDS(combined, file = paste0(results_path, "combined_clustree.rds"))

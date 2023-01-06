#'velocity-extract from Seurat
#'Data: 2022/6/5
#'Author: Han Zhifa
#'Reference: https://github.com/basilkhuder/Seurat-to-RNA-Velocity

####packages and workDir####
try(library(Seurat, lib.loc = "/home/hanzf/.conda/envs/scs10x/lib/R/library"), 
    silent = T)
require(Seurat, lib.loc = "/opt/R/4.1.0/lib/R/library")
#getwd()
setwd("/ssd/hanzf/research/scs10x")
####get data####
# get results outputing path
results_path <- paste0(commandArgs(T)[1], "extracted_data/")
dir.create(results_path)
# read combined
combined <- readRDS(commandArgs(T)[2])
#### Extracting Meta-data ####
# Filtered Cell Ids
write.csv(Cells(combined), 
          file = paste0(results_path, "cellID_obs.csv"), 
          row.names = FALSE)
# UMAP or TSNE coordinates
write.csv(Embeddings(combined, reduction = "umap"), 
          file = paste0(results_path, "cell_embeddings.csv"))
# Clusters (Optional)
write.csv(data.frame(CellID = Cells(combined), clusters = combined@meta.data$seurat_clusters), 
          file = paste0(results_path, "clusters.csv"))
# celltypes (Optional)
write.csv(data.frame(CellID = Cells(combined), celltypes = combined@meta.data$celltypes), 
          file = paste0(results_path, "celltypes.csv"))

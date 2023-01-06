#' Title: make loom file
#' Author: Han Zhifa
#' Date: 2021.11.19
#' Reference1: https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Setup.html#c-from-other-r-objects-e.g.-seurat-singlecellexperiment
#' Reference2: https://htmlpreview.github.io/?https://github.com/aertslab/SCopeLoomR/blob/master/vignettes/SCopeLoomR_Seurat_tutorial.nb.html
#' 
try(library(Seurat, lib.loc = "/home/hanzf/.conda/envs/scs10x/lib/R/library"), 
    silent = T)
require(Seurat, lib.loc = "/opt/R/4.1.0/lib/R/library")
library(SCopeLoomR)
library(SCENIC)
setwd("/ssd/hanzf/research/scs10x/")
dir.create(paste0(commandArgs(T)[2], "SCENIC/"))

combined <- readRDS(commandArgs(T)[1])

# Default
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltypes"

## loom
exprMat <- combined@assays$RNA@counts
print( colnames(combined@meta.data) )
cellInfo <- combined@meta.data[, c("orig.ident", "Phase", 
                                   "experiment", "celltypes")]

loom <- build_loom(paste0(commandArgs(T)[2], "SCENIC/", "combined_filtered.loom"),
                   dgem=exprMat, title = "combined", 
                   default.embedding = combined@reductions$umap@cell.embeddings, 
                   default.embedding.name = "umap.rna")

## add clustering
loom <- open_loom(paste0(commandArgs(T)[2], "SCENIC/", "combined_filtered.loom"), mode = "r+")
add_clustering(loom = loom,
               group = "Custom clustering",
               name = "celltypes",
               clusters = combined$celltypes, 
               #annotation = combined$celltypes, 
               #overwrite.default = T,
               is.default = T )

#add_seurat_clustering(loom = loom,
#                      seurat = combined,
#                      seurat.assay = "integrated",
#                      seurat.clustering.prefix = "integrated_snn_res.",
#                      default.clustering.resolution = "res.0.3",
#                      )

## add meta data
while (F) {
  # discrete attritte
  Attri_discrete <- c("orig.ident", "experiment", "doublet", "seurat_clusters", "Phase", 
                      "celltypes")
  for (ad in Attri_discrete) {
    add_col_attr(loom = loom, key = paste0(ad, "2"), 
                 value = combined@meta.data[, ad], as.annotation = TRUE)
  }
  # continues attritte
  Attri_continues <- c("sum")
  for (ac in Attri_continues) {
    add_col_attr(loom = loom, key = ac,
                 value = combined@meta.data[, ac], as.metric = TRUE)
  }
}

loom <- add_cell_annotation(loom, cellInfo)

# Always remember to close loom files when done
close_loom(loom)




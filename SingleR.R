#'SingleR
#'Author: Han Zhifa
#'date: Jan 8th 2020
#*####**************************************####*#
#             packages and workDir               #
######======================================######
.libPaths()
#BiocManager::install(version = "3.10", lib = "/usr/lib64/R/library")
#BiocManager::install("SingleR")
#BiocManager::install("SingleCellExperiment")
#BiocManager::install("scater")

library(SingleR)
library(scater)
try(library(Seurat, lib.loc = "/opt/R/4.1.0/lib/R/library"), 
    silent = T)
require(Seurat, lib.loc = "/home/hanzf/.conda/envs/scs10x/lib/R/library")
library(SingleCellExperiment)
library(celldex)

setwd("/ssd/hanzf/research/scs10x")

#####

#*####**************************************####*#
#                    data                        #
######======================================######
## 
results_path <- commandArgs(T)[1]

## read combined
combined <- readRDS(commandArgs(T)[2])

## create marker path
singler_path <- paste0(results_path, "SingleR/")
dir.create(singler_path)

#####

#*####**************************************####*#
#                   SingleR                      #
######======================================######
sce <- as.SingleCellExperiment(combined)

#sce <- logNormCounts(sce)

## reference database
hpca.se <- HumanPrimaryCellAtlasData()
#IG <- ImmGenData()
#MR <- MouseRNAseqData()
DICE <- DatabaseImmuneCellExpressionData()
NHD <- NovershternHematopoieticData()
MID <- MonacoImmuneData()

ref <- list(HPCA = hpca.se, DICE = DICE, 
            NHD = NHD, MID = MID)

## annotate the celltypes
for(i in 1:length(ref)){
  celltypes <- SingleR(test = sce, ref = ref[[i]], labels = ref[[i]]$label.main)
  addmeta <- data.frame(SingleR.labels = celltypes@listData$labels, 
                        row.names = celltypes@rownames)
  
  combined <- AddMetaData(combined, addmeta)
  
  png(paste0(singler_path, names(ref)[i], "_UMAP.png"), 
      width = 480 * 2, height = 480 * 2)
  print(DimPlot(combined, reduction = "umap", label = T, group.by = "SingleR.labels"))
  dev.off()
  
  jpeg(paste0(singler_path, names(ref)[i], "_violin.jpeg"), quality = 100, width = 1600, 
       height = ceiling(length(unique(celltypes$labels)) / 3) * 600)
  print(plotScoreDistribution(celltypes, show = "delta.med", ncol = 3, show.nmads = 3, 
                              size = 0.2))
  dev.off()
  
  pdf(paste0(singler_path, names(ref)[i], "_heatmap.pdf"))
  plotScoreHeatmap(celltypes)
  dev.off()
  
  write.csv(addmeta, 
            file = paste0(singler_path, names(ref)[i], ".csv"))
}

if(F){
  
  celltypes <- SingleR(test = sce, ref = MR, labels = MR$label.main)
  addmeta <- data.frame(SingleR.labels = celltypes@listData$labels, 
                        row.names = celltypes@rownames)
  
  combined <- AddMetaData(combined, addmeta)
  SingleR.labels
  
  png(paste0(singler_path, i, "_UMAP.png"), 
      width = 480 * 2, height = 480 * 2)
  DimPlot(combined, reduction = "umap", label = T, group.by = "SingleR.labels")
  dev.off()
  
  pdf()
  plotScoreHeatmap(celltypes)
  dev.off()
  
  pdf()
  plotScoreHeatmap(celltypes, 
                   annotation_col=as.data.frame(colData(sceG)[,"donor",drop=FALSE]))
  dev.off()
  
  pdf()
  plotScoreDistribution(celltypes, show = "delta.med", ncol = 3, show.nmads = 3)
  dev.off()
  
  all.markers <- metadata(celltypes)$de.genes
  sce$labels <- celltypes$labels
  
  # Beta cell-related markers
  plotHeatmap(sce, order_columns_by="labels",
              features=unique(unlist(all.markers$beta)))
  
}

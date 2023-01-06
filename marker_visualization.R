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
#install.packages("pheatmap")
#install.packages("foreach")
try(library(Seurat, lib.loc = "/home/hanzf/.conda/envs/scs10x/lib/R/library"), 
    silent = T)
require(Seurat, lib.loc = "/opt/R/4.1.0/lib/R/library")
#require(future)
require(foreach)
#require(cowplot)

#getwd()
setwd("/ssd/hanzf/research/scs10x")

## set parallelization
#cat("set the workers used in visualization (according cores): ")
#workers <- as.numeric(readLines("stdin", n = 1))
#plan("multiprocess", workers = workers)
# 
#plan()

#####

#*####**************************************####*#
#                    get data                    #
######======================================######
## 
results_path <- commandArgs(T)[1]

## read combined
combined <- readRDS(commandArgs(T)[2])

## read marker filenames
ref_marker_file <- dir(commandArgs(T)[3])

## read all_markers
all_markers <- read.csv(commandArgs(T)[4])

## create marker path
marker_path <- paste0(results_path, "ref_markers_VSL/")
dir.create(marker_path)

#####

#*####**************************************####*#
#         all markers and cluster markers        #
######======================================######
##
DefaultAssay(combined) <- "RNA"

cat(paste0("Do you want to plot RidgePlot (take a long time) (T or F):"))
wheRP <- as.logical(readLines("stdin", n=1))
## identify the ref markers
foreach(f = ref_marker_file) %do% {
  ## output cluster markers
  ref_markers <- read.table(paste0("results/ref_markers/", f), header = T, sep = "\t")
  cluster_markers <- merge(ref_markers, all_markers, by = "gene")
  cluster_markers <- cluster_markers[order(cluster_markers$cluster), ]
  write.table(cluster_markers, file = paste0(marker_path, "ref_", f), sep = "\t", 
              row.names = F)
  
  ## RidgePlot
  if(wheRP & length(unique(cluster_markers$gene)) < 190){
    pdf(paste0(marker_path, "RidegPlot_ref_", f, ".pdf"), 
        height = length(unique(combined$seurat_clusters)) * 0.09 * 
          length(unique(cluster_markers$gene)) + 0.5, 
        width = 19)
    print(RidgePlot(combined, features = as.character(unique(cluster_markers$gene))))
    dev.off()
  
  
  
  ## DotPlot
  if(length(unique(cluster_markers$gene)) < 190){
    pdf(paste0(marker_path, "DotPlot_ref_", f, ".pdf"), 
        width = length(unique(cluster_markers$gene))* 0.17 + 5, 
        height = length(unique(combined$seurat_clusters)) * 0.27 + 3)
    print(DotPlot(combined, features = rev(as.character(unique(cluster_markers$gene))), 
                  group.by = "seurat_clusters") + RotatedAxis())
    dev.off()
  }
  
  ## VlnPlot
  if(length(unique(cluster_markers$gene)) < 190){
    jpeg(paste0(marker_path, "VlnPlot_ref_", f, ".jpeg"), quality = 100, 
         width = length(unique(combined$seurat_clusters)) * 45 + 200, 
         height = ceiling(length(unique(cluster_markers$gene)) / 3) * 450)
    print(VlnPlot(object = combined, ncol = 3, 
                  features = as.character(unique(cluster_markers$gene))))
    dev.off()
  }
  }
}

foreach(f = ref_marker_file) %do% {
  ## output cluster markers
  ref_markers <- read.table(paste0("results/ref_markers/", f), header = T, sep = "\t")
  cluster_markers <- merge(ref_markers, data.frame(gene = rownames(combined)), 
                           by = "gene")
  print(colnames(cluster_markers))
  
  ## get marker gene used for ploting
  cluster_markers <- cluster_markers[order(cluster_markers$celltype), ]
  
  ## sub !\w character in celltype names
  cluster_markers$celltype <- gsub(pattern = "[^\\s\\w]", replacement = " ", 
                                   cluster_markers$celltype, perl = T)
  
  ## creat dir and plot
  dir.create(paste0(marker_path, "DotPlot_", f))
  dir.create(paste0(marker_path, "RidegPlot_", f))
  dir.create(paste0(marker_path, "VlnPlot_", f))
  with(cluster_markers, foreach(i = unique(celltype)) %do% {
    num <- i == celltype
    
    ## DotPlot
    if(length(unique(gene[num])) < 199){
      pdf(paste0(marker_path, "DotPlot_", f, "/", i, ".pdf"), 
          width = length(unique(gene[num]))* 0.17 + 5, 
          height = length(unique(combined$seurat_clusters)) * 0.27 + 3)
      print(DotPlot(combined, features = rev(as.character(unique(gene[num]))), 
                    group.by = "seurat_clusters") + RotatedAxis())
      dev.off()
    }else{
      pdf(paste0(marker_path, "DotPlot_", f, "/", i, ".pdf"), 
          width = length(unique(gene[num])[1:199])* 0.17 + 5, 
          height = length(unique(combined$seurat_clusters)) * 0.27 + 3)
      print(DotPlot(combined, features = rev(as.character(unique(gene[num][1:199]))), 
                    group.by = "seurat_clusters") + RotatedAxis())
      dev.off()
    }
    
    ## RidgePlot
    if(wheRP){
      if(length(unique(gene[num])) < 199){
        pdf(paste0(marker_path, "RidegPlot_", f, "/", i, ".pdf"), 
            height = length(unique(combined$seurat_clusters)) * 0.09 * 
              length(unique(gene[num])) + 0.5, 
            width = 19)
        print(RidgePlot(combined, features = as.character(unique(gene[num]))))
        dev.off()
      }else{
        pdf(paste0(marker_path, "RidegPlot_", f, "/", i, ".pdf"), 
            height = length(unique(combined$seurat_clusters)) * 0.09 * 
              length(unique(gene[num][1:199])) + 0.5, 
            width = 19)
        print(RidgePlot(combined, features = as.character(unique(gene[num][1:199]))))
        dev.off()
      }
    
    
    ## VlnPlot
    if(length(unique(gene[num])) < 199){
      jpeg(paste0(marker_path, "VlnPlot_", f, "/", i, ".jpeg"), quality = 100, 
           width = length(unique(combined$seurat_clusters)) * 45 + 200, 
           height = ceiling(length(unique(gene[num])) / 3) * 450)
      print(VlnPlot(object = combined, 
                    features = as.character(unique(gene[num])), 
                    ncol = 3))
      dev.off()
    }else{
      jpeg(paste0(marker_path, "VlnPlot_", f, "/", i, ".jpeg"), quality = 100, 
           width = length(unique(combined$seurat_clusters)) * 45 + 200, 
           height = ceiling(length(unique(gene[num][1:199])) / 3) * 450)
      print(VlnPlot(object = combined, 
                    features = as.character(unique(gene[num][1:199])), 
                    ncol = 3))
      dev.off()
    }
    }
  })
}

#####

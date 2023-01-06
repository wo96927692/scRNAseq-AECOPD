#'10x data: intergrate different datasets
#'Data: 2019/09/26
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
try(library(Seurat, lib.loc = "/opt/R/4.1.0/lib/R/library"), 
    silent = T)
require(Seurat, lib.loc = "/home/hanzf/.conda/envs/scs10x/lib/R/library")

require(pheatmap)

#getwd()
setwd("/ssd/hanzf/research/scs10x")

#####

#*####**************************************####*#
#                    get data                    #
######======================================######
## input the results path for outputing is necessary
results_path <- commandArgs(T)[1]
dir.create(results_path)
PC_path <- paste0(results_path, "PCassess/")
dir.create(PC_path)

## read combined
combined <- readRDS(commandArgs(T)[2])

## read split_factor
split_factor <- as.character(read.csv(commandArgs(T)[3], header = T, quote = "\"")$split_factor)
if(split_factor == "NULL") split_factor <- NULL

#####

#*####**************************************####*#
#              dimensional reduction             #
######======================================######
## PCA
combined <- RunPCA(combined, npcs = 50)
# Examine and visualize PCA results a few different ways
png(paste0(PC_path, "PCA.png"), 
    width = 480*2, height = 480*2)
DimPlot(object = combined, reduction = "pca", group.by = split_factor)
dev.off()

## tSNE
combined <- RunTSNE(combined, dims.use = 1:10, do.fast = TRUE, 
                    check_duplicates = FALSE)
# Examine and visualize PCA results a few different ways
png(paste0(PC_path, "tSNE.png"), 
    width = 480*2, height = 480*2)
DimPlot(object = combined, reduction = "tsne", group.by = split_factor)
dev.off()


## UMAP
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
# Examine and visualize PCA results a few different ways
png(paste0(PC_path, "UMAP.png"), 
    width = 480*2, height = 480*2)
DimPlot(object = combined, reduction = "umap", group.by = split_factor)
dev.off()

#####

#*####**************************************####*#
#                 assessing PCs                  #
######======================================######
## assessing PCs
## JackStraw
combined <- JackStraw(object = combined, reduction = "pca", dims = 50, 
                      num.replicate = 100,  prop.freq = 0.1, verbose = 1)
combined <- ScoreJackStraw(object = combined, dims = 1:50, reduction = "pca")
# plot JackStraw
pdf(paste0(PC_path, "JackStraw.pdf"), width = 50/20 * 1.5 + 10)
JackStrawPlot(object = combined, dims = 1:50, reduction = "pca")
dev.off()

## plot the sd of every PCs
pdf(paste0(PC_path, "Elbow.pdf"))
ElbowPlot(object = combined)
dev.off()

## plot the heatmap of PCs vs Genes
seq <- c(seq(1, 15, 3), 16:24, seq(25, 50, 5))
for(i in seq){
  pdf(paste0(PC_path, "PCA_Heatmap_PC", i, ".pdf"))
  print(DimHeatmap(object = combined, reduction = "pca", nfeatures = 50, 
                   cells = 200, balanced = TRUE, dims = i))
  dev.off()
  cat(i)
}

#####

#*####**************************************####*#
#                     get PCs                    #
######======================================######
# read PCs
cat("select the PCs used in download analysis (sperate by ',' (eg: 1,2,3) or use ':' operation (eg: 1:3)): \n")
PCs <- readLines("stdin", n = 1)
if(grepl(":", PCs)){
  PCs <- eval(parse(text = PCs))
}else{
  PCs <- as.numeric(unlist(strsplit(PCs, ",")))
}
cat("The used PCs contain: "); cat(PCs)

#####

## save data
write.csv(data.frame(PCs = PCs), file = paste0(results_path, "PCs.csv"), 
          row.names = F)
saveRDS(combined, file = paste0(results_path, "combined_PCassessing.rds"))

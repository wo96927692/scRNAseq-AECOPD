#' Title: Comparison analysis (CellChat)
#' Author: Han Zhifa
#' Date: 2021.10.08
#' 
#*####**************************************####*#
#       packages and workDir and data            #
######======================================######
.libPaths("/home/hanzf/.conda/envs/scs10x/lib/R/library")
library(CellChat)
library(patchwork)
library(ComplexHeatmap)

setwd("/ssd/hanzf/research/scs10x/")

results_path <- commandArgs(T)[1]
CellChat_path <- paste0(results_path, "CellChat_COMP/")
dir.create(CellChat_path)

cat("please input the absolute or relative filepath for cellchat object 1: \n")
file_cellchat_1 <- readLines("stdin", n = 1)
cat("please input the absolute or relative filepath for cellchat object 2: \n")
file_cellchat_2 <- readLines("stdin", n = 1)

cellchat_CT <- readRDS(file_cellchat_1)
cellchat_CS <- readRDS(file_cellchat_2)

object.list <- list(CT = cellchat_CT, CS = cellchat_CS)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#### Part I: comparison ####
#barplot
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), 
                           measure = "weight")
pdf(file = paste0(CellChat_path, "Barplot comparison.pdf"), 
    height = 4)
gg1 + gg2
dev.off()
#circle plot
pdf(file = paste0(CellChat_path, "Circleplot comparison.pdf"))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()
#Heatmap
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
pdf(file = paste0(CellChat_path, "Heatmap comparison.pdf"), 
    width = 14)
gg1 + gg2
dev.off()

##Circle plot for each dateset
weight.max <- getMaxWeight(object.list, 
                           attribute = c("idents","count"))
pdf(file = paste0(CellChat_path, "Circleplot each.pdf"), 
    width = 14)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, 
                   weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], 
                   edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", 
                                       names(object.list)[i])
                   )
}
dev.off()
# Scatter for each dateset
num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)
  })
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- 
    netAnalysis_signalingRole_scatter(object.list[[i]], 
                                      title = names(object.list)[i], 
                                      weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf(file = paste0(CellChat_path, "Scatter each.pdf"), 
    width = 14)
patchwork::wrap_plots(plots = gg)
dev.off()

#### Part II: signaling pathways ####
## functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, 
                                         type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
pdf(file = paste0(CellChat_path, "estimationNumCluster_functional.pdf"))
cellchat <- netClustering(cellchat, type = "functional")
dev.off()
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
pdf(file = paste0(CellChat_path, "functional similarity.pdf"))
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
dev.off()
# each cluster
pdf(file = paste0(CellChat_path, 
                  "functional similarity each cluster.pdf"))
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", 
                                  nCol = 2)
#> 2D visualization of signaling networks from datasets 1 2
dev.off()
# distance
pdf(file = paste0(CellChat_path, "functional distance.pdf"))
rankSimilarity(cellchat, type = "functional")
#> Compute the distance of signaling networks between datasets 1 2
dev.off

## structural similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
pdf(file = paste0(CellChat_path, "estimationNumCluster_structural.pdf"))
cellchat <- netClustering(cellchat, type = "structural")
dev.off()
#> Classification learning of the signaling networks for datasets 1 2
pdf(file = paste0(CellChat_path, "structural similarity.pdf"))
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
dev.off()
# each cluster
pdf(file = paste0(CellChat_path, 
                  "structural similarity each cluster.pdf"))
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
#> 2D visualization of signaling networks from datasets 1 2
dev.off()
# distance
pdf(file = paste0(CellChat_path, "structural distance.pdf"))
rankSimilarity(cellchat, type = "structural")
#> Compute the distance of signaling networks between datasets 1 2
dev.off

## information flow
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
pdf(file = paste0(CellChat_path, "information flow.pdf"), width = 14)
gg1 + gg2
dev.off()

## outgoing (or incoming) signaling vs each cell population
# outgoing
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, 
                       object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], 
                                        pattern = "outgoing", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i], 
                                        width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], 
                                        pattern = "outgoing", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i+1], 
                                        width = 5, height = 6)
pdf(file = paste0(CellChat_path, "outgoing signals and celltypes.pdf"), 
    width = 12)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
# incombing
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], 
                                        pattern = "incoming", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i], 
                                        width = 5, height = 6, 
                                        color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], 
                                        pattern = "incoming", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i+1], 
                                        width = 5, height = 6, 
                                        color.heatmap = "GnBu")
pdf(file = paste0(CellChat_path, "incoming signals and celltypes.pdf"), 
    width = 12)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
#
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], 
                                        pattern = "all", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i], 
                                        width = 5, height = 6, 
                                        color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], 
                                        pattern = "all", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i+1], 
                                        width = 5, height = 6, 
                                        color.heatmap = "OrRd")
pdf(file = paste0(CellChat_path, 
                  "all signals and celltypes.pdf"), width = 12)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

## save data
saveRDS(cellchat, file = paste0(CellChat_path, 
                                "cellchat_comparison.rds"))

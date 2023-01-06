#' Title: Customization (CellChat)
#' Author: Han Zhifa
#' Date: 2021.10.08
#' 
#*####**************************************####*#
#       packages and workDir and data            #
######======================================######
.libPaths("/home/hanzf/.conda/envs/scs10x/lib/R/library")
library(CellChat)
library(ggplot2, 
        lib.loc = c("/usr/lib64/R/library", "/usr/share/R/library") )
library(ggalluvial)
library(NMF)
library(Seurat, lib.loc = c("/opt/R/4.1.0/lib/R/library"))
library(patchwork)

setwd("/ssd/hanzf/research/scs10x/")


#*####**************************************####*#
#            Comparison analysis                 #
######======================================######
# read data
cellchat <- readRDS("results/cluster/CellChat_COMP/cellchat_comparison.rds")

cellchat_CT <- readRDS("results/cluster/CellChat/experiment/Stable COPD/ALL/cellchat_ALL.rds")
cellchat_CS <- readRDS("results/cluster/CellChat/experiment/AECOPD/ALL/cellchat_ALL.rds")
object.list <- list(CT = cellchat_CT, CS = cellchat_CS)


## Part I: Differential number of interactions or interaction strength among different cell types ####
print( levels (cellchat@idents$joint) )
# 
cell1 <- "CD14+Mono"; cell2 <- "NK"; cell3 <- "NKG2C+NK"
group.cellType <- c(rep(cell1, 4), rep(cell2, 4), rep(cell3, 4))
group.cellType <- factor(group.cellType, levels = c(cell1, cell2, cell3))
object.list <- lapply(object.list, function(x) {
  mergeInteractions(x, group.cellType)
  })
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
# plot
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, 
                          measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, 
                          measure = "weight.merged", label.edge = T)

## sources and targets in specific celltype
print( levels (cellchat@idents$joint) )
celltype <- "NKG2C+NK" # choose one celltype for visualization
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, 
                                            idents.use = celltype, 
                                            signaling.exclude = "IGTB2")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1
# plot
patchwork::wrap_plots(plots = list(gg1))

#### Part III: upgulated and down-regulated signaling ####
print( levels (cellchat@idents$joint) )
sources_celltypes <- c(3, 11) # rank number
targets_celltypes <- c(4) # rank number
netVisual_bubble(cellchat, sources.use = sources_celltypes, 
                 targets.use = targets_celltypes, 
                 comparison = c(1, 2), 
                 angle.x = 45)
#> Comparing communications on a merged object

gg1 <- netVisual_bubble(cellchat, sources.use = sources_celltypes, 
                        targets.use = targets_celltypes,  
                        comparison = c(1, 2), max.dataset = 2, 
                        title.name = "Increased signaling in CS", 
                        angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = sources_celltypes, 
                        targets.use = targets_celltypes,  
                        comparison = c(1, 2), max.dataset = 1, 
                        title.name = "Decreased signaling in CS", 
                        angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

#### Expression analysis
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "CS"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, 
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.1, 
                                       thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "CS",
                              ligand.logFC = 0.2, 
                              receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, 
                                datasets = "CT",
                                ligand.logFC = -0.1, 
                                receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, 
                        sources.use = 4, targets.use = c(5:11), 
                        comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,
                        title.name = 
                          paste0("Up-regulated signaling in ", 
                                 names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, 
                        sources.use = 4, targets.use = c(5:11), 
                        comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,title.name = 
                          paste0("Down-regulated signaling in ", 
                                 names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], 
                     sources.use = sources_celltypes, 
                     targets.use = targets_celltypes, 
                     slot.name = 'net', 
                     net = net.up, lab.cex = 0.8, small.gap = 3.5, 
                     title.name = 
                       paste0("Up-regulated signaling in ", 
                              names(object.list)[2]))
#> Note: The first link end is drawn out of sector 'MIF'.
netVisual_chord_gene(object.list[[1]], 
                     sources.use = sources_celltypes, 
                     targets.use = targets_celltypes, 
                     slot.name = 'net', 
                     net = net.down, lab.cex = 0.8, small.gap = 3.5, 
                     title.name = 
                       paste0("Down-regulated signaling in ", 
                                         names(object.list)[2]))


#### Part IV: specific pathway ####
pathways.show <- c("MHC-I") 

## circle plot
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), 
                           attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, 
                      layout = "circle", edge.weight.max = 
                        weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, 
                                             names(object.list)[i]))
}

## heatmap
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], 
                               signaling = pathways.show, 
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, 
                                                  "signaling ",
                                                  names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

## Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, 
                      layout = "chord", 
                      signaling.name = paste(pathways.show, 
                                             names(object.list)[i]))
}
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

## Chord diagram 2
group.cellType <- c(rep(cell1, 4), rep(cell2, 4), rep(cell3, 4)) # grouping cell clusters into fibroblast, DC and TC cells names(group.cellType) <- levels(object.list[[1]]@idents)
names(group.cellType) <- levels(object.list[[1]]@idents)

par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, 
                                           " signaling network - ", 
                                           names(object.list)[i]))
}
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], 
                       sources.use = sources_celltypes, 
                       targets.use = targets_celltypes, 
                       #lab.cex = 0.5, 
                       title.name = 
                         paste0("Signaling in ", 
                                names(object.list)[i]))
}

# compare all the interactions sending from fibroblast to inflamatory immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], 
                       sources.use = sources_celltypes, 
                       targets.use = targets_celltypes,  
                       title.name = 
                         paste0("Signaling in ", 
                                names(object.list)[i]), 
                       legend.pos.x = 10)
}

# show all the significant signaling pathways from fibroblast to immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], 
                       sources.use = sources_celltypes, 
                       targets.use = targets_celltypes,
                       slot.name = "netP", 
                       title.name = 
                         paste0("Signaling pathways in ", 
                                names(object.list)[i]), 
                       legend.pos.x = 10)
}
#> Note: The second link end is drawn out of sector ' '.
#> Note: The first link end is drawn out of sector 'MIF'.
#> Note: The second link end is drawn out of sector ' '.
#> Note: The first link end is drawn out of sector 'CXCL '.

#### Part V: Genes ####
cellchat@meta$datasets = factor(cellchat@meta$datasets, 
                                levels = c("CT", "CS")) # set factor level
plotGeneExpression(cellchat, signaling = pathways.show, 
                   split.by = "datasets", 
                   colors.ggplot = T)


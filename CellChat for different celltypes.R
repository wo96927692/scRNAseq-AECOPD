#' Title: cell cell communication (CellChat) for different NK
#' Author: Han Zhifa
#' Date: 2021.09.25
#' 
#*####**************************************####*#
#       packages and workDir and data            #
######======================================######
library(CellChat, lib.loc = "/opt/R/4.1.0/lib/R/library")
library(ggalluvial)
try(library(Seurat, lib.loc = "/home/hanzf/.conda/envs/scs10x/lib/R/library"), 
    silent = T)
require(Seurat, lib.loc = "/opt/R/4.1.0/lib/R/library")
library(NMF)
library(patchwork)
library(ComplexHeatmap)


setwd("/ssd/hanzf/research/scs10x/")
#### I: ####
combined_all <- readRDS("results/cluster/combined_celltypes.rds")
Idents(combined_all) <- "celltypes"

species <- "human"

results_path <- "results/cluster/"
CellChat_path_0 <- paste0(results_path, "CellChat_all/")
dir.create(CellChat_path_0)


CellChat_path_1 <- paste0(CellChat_path_0, "celltypes/")
#CellChat_path_1 <- CellChat_path_0
dir.create(CellChat_path_1)


if(species == "human"){
  ## Set the ligand-receptor interaction database
  CellChatDB <- CellChatDB.human
  pdf(file = paste0(CellChat_path_0, "DB.pdf"))
  print(showDatabaseCategory(CellChatDB))
  dev.off()
  # Show the structure of the database
  print(dplyr::glimpse(CellChatDB$interaction))
}else{
  ## Set the ligand-receptor interaction database
  CellChatDB <- CellChatDB.mouse
  pdf(file = paste0(CellChat_path_0, "DB.pdf"))
  print(showDatabaseCategory(CellChatDB))
  dev.off()
  # Show the structure of the database
  print(dplyr::glimpse(CellChatDB$interaction))
}

e = "NK" # celltype you want to exclude
old_celltype <- "NKG2C+NK"; new_celltype <- "NK"
#CellChat_path <- CellChat_path_0
CellChat_path <- paste0(paste0(CellChat_path_1, "no-", e, "/"))
dir.create(CellChat_path, recursive = T)

##### creat cellchat object
combined <- combined_all
combined <- combined_all[, combined_all$celltypes != e ]
celltype_rank = which(levels(combined$celltypes) == old_celltype)
new_celltypes <- levels(combined$celltypes)
new_celltypes[celltype_rank] <- new_celltype
levels(combined$celltypes) <- new_celltypes
Idents(combined) <- "celltypes"
data.input <- GetAssayData(combined, assay = "RNA", 
                           slot = "data") # normalized data matrix
labels <- Idents(combined)
meta <- data.frame(labels = labels, 
                   row.names = names(labels)) # create a dataframe of the cell labels

print(unique(meta$labels)) # check the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
print(levels(cellchat@idents)) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

##### cell-cell communication analysis ####
i ="ALL"
subCellChat_path <- paste0(CellChat_path, i, "/")
dir.create(subCellChat_path)
subcellchat <- cellchat

if(i == "ALL"){
  CellChatDB.use <- CellChatDB
}else{
  CellChatDB.use <- subsetDB(CellChatDB, search = i)# "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
}
subcellchat@DB <- CellChatDB.use # set the used database in the object

subcellchat <- subsetData(subcellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
subcellchat <- identifyOverExpressedGenes(subcellchat)
subcellchat <- identifyOverExpressedInteractions(subcellchat)
subcellchat <- projectData(subcellchat, PPI.human)
## Inference of cell-cell communication network
subcellchat <- computeCommunProb(subcellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
subcellchat <- filterCommunication(subcellchat, min.cells = 10)
subcellchat <- computeCommunProbPathway(subcellchat)
subcellchat <- aggregateNet(subcellchat)

## visualization
# All circle
groupSize <- as.numeric(table(subcellchat@idents))
pdf(file = paste0(subCellChat_path, "All_circle.pdf"))
par(mfrow = c(1,2), xpd=TRUE)
print(netVisual_circle(subcellchat@net$count, vertex.weight = groupSize, 
                       weight.scale = T, label.edge= F, 
                       title.name = "Number of interactions"))
print(netVisual_circle(subcellchat@net$weight, vertex.weight = groupSize, 
                       weight.scale = T, label.edge= F, 
                       title.name = "Interaction weights/strength"))
dev.off()
# Sub circle
pdf(file = paste0(subCellChat_path, "SUB_circle.pdf"))
mat <- subcellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (j in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), 
                 dimnames = dimnames(mat))
  mat2[j, ] <- mat[j, ]
  print(netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                         edge.weight.max = max(mat), 
                         title.name = rownames(mat)[j]))
}
dev.off()
# Access all the signaling pathways showing significant communications
pathways.show.all <- subcellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(subcellchat@idents)
vertex.receiver = seq(1,4)
DFWD <- getwd()
dir.create(paste0(subCellChat_path, "signal/"))
setwd(paste0(subCellChat_path, "signal/"))
for (k in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(subcellchat, signaling = pathways.show.all[k], 
            vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(subcellchat, signaling = pathways.show.all[k])
  ggsave(filename=paste0(pathways.show.all[k], 
                         "_L-R_contribution.png"), plot=gg, width = 6, 
         height = 4, units = 'in', dpi = 300)
}
setwd(DFWD)
# All bubble
pdf(file = paste0(subCellChat_path, "All_bubble.pdf"), 
    height = 12, width = 25)
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
print(netVisual_bubble(subcellchat, 
                       sources.use = 1:length(levels(subcellchat@idents)), 
                       targets.use = 1:length(levels(subcellchat@idents)), 
                       remove.isolate = FALSE))
#> Comparing communications on a single object
dev.off()

# All Chord
pdf(file = paste0(subCellChat_path, "All_chord.pdf"), 
    width = 10)
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
print(netVisual_chord_gene(subcellchat, 
                           sources.use = 1:length(levels(subcellchat@idents)), 
                           targets.use = 1:length(levels(subcellchat@idents)), 
                           lab.cex = 0.5,legend.pos.y = 50))
dev.off()

#### Systems analysis ####
# Compute the network centrality scores
subcellchat <- netAnalysis_computeCentrality(subcellchat, 
                                             slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
#
pdf(file = paste0(subCellChat_path, "Domaint senders and receivers.pdf"))
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
print(netAnalysis_signalingRole_scatter(subcellchat))
dev.off()
#
pdf(file = paste0(subCellChat_path, 
                  "Celltypes contributions of senders and receivers.pdf"), 
    width = 15)
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(subcellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(subcellchat, pattern = "incoming")
print(ht1 + ht2)
dev.off()

## outgoing
pdf(file = paste0(subCellChat_path, 
                  "outgoing nPattern decision.pdf"), 
    height = 4)
print(selectK(subcellchat, pattern = "outgoing"))
dev.off()
#
# the number of patterns of outgoing
nPatterns <- 4
#
pdf(file = paste0(subCellChat_path, "outgoing Communication Patterns.pdf"), 
    height = 4, width = 10)
cellchat <- identifyCommunicationPatterns(subcellchat, 
                                          pattern = "outgoing", 
                                          k = nPatterns)
dev.off()
#
pdf(file = paste0(subCellChat_path, "outgoing river.pdf"))
# river plot
try(print(netAnalysis_river(subcellchat, pattern = "outgoing")), silent = T)
#> Please make sure you have load `library(ggalluvial)` when running this function
dev.off()
# dot plot
pdf(file = paste0(subCellChat_path, "outgoing dot.pdf"))
try(print(netAnalysis_dot(cellchat, pattern = "outgoing")), silent = T)
dev.off()

## incoming
pdf(file = paste0(subCellChat_path, 
                  "incoming nPattern decision.pdf"), 
    height = 4)
print(selectK(subcellchat, pattern = "incoming"))
dev.off()
#
nPatterns <- 5
#
pdf(file = paste0(subCellChat_path, "incoming Communication Patterns.pdf"), 
    height = 8, width = 10)
subcellchat <- identifyCommunicationPatterns(subcellchat, 
                                             pattern = "incoming", 
                                             k = nPatterns)
dev.off()
#
pdf(file = paste0(subCellChat_path, "incoming river.pdf"), 
    height = 6)
# river plot
try(print(netAnalysis_river(subcellchat, pattern = "incoming")), silent = T)
#> Please make sure you have load `library(ggalluvial)` when running this function
dev.off()
# dot plot
pdf(file = paste0(subCellChat_path, "incoming dot.pdf"))
try(print(netAnalysis_dot(cellchat, pattern = "incoming")), silent = T)
dev.off()

setwd(subCellChat_path)
subcellchat <- computeNetSimilarity(subcellchat, type = "functional")
subcellchat <- netEmbedding(subcellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
subcellchat <- netClustering(subcellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
setwd(DFWD)
# Visualization in 2D-space
pdf(file = paste0(subCellChat_path, "netVisual functional.pdf"))
print(netVisual_embedding(subcellchat, type = "functional", 
                          label.size = 3.5))
dev.off()
#
pdf(file = paste0(subCellChat_path, "netVisual functional zoomin.pdf"))
netVisual_embeddingZoomIn(subcellchat, type = "functional", nCol = 2)
dev.off()

## 
setwd(subCellChat_path)
subcellchat <- computeNetSimilarity(subcellchat, type = "functional")
subcellchat <- computeNetSimilarity(subcellchat, type = "structural")
subcellchat <- netEmbedding(subcellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
subcellchat <- netClustering(subcellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
setwd(DFWD)
# Visualization in 2D-space
pdf(file = paste0(subCellChat_path, "netVisual structure.pdf"))
netVisual_embedding(subcellchat, type = "structural", label.size = 3.5)
dev.off()
#
pdf(file = paste0(subCellChat_path, "netVisual structure zoomin.pdf"))
netVisual_embeddingZoomIn(subcellchat, type = "structural", nCol = 2)
dev.off()

## save data
saveRDS(subcellchat, 
        file = paste0(subCellChat_path, "cellchat_", i, ".rds"))


#### II: Comparison #### 
results_path <- "results/cluster/"
CellChat_path <- paste0(CellChat_path_1, "CellChat_COMP/")
dir.create(CellChat_path)

# read the file for cellchat object 1)
file_cellchat_1 <- paste0(CellChat_path_1, "no-NKG2C+NK/ALL/cellchat_ALL.rds")
# read the file for cellchat object 2)
file_cellchat_2 <- paste0(CellChat_path_1, "no-NK/ALL/cellchat_ALL.rds")

cellchat_CT <- readRDS(file_cellchat_1)
cellchat_CS <- readRDS(file_cellchat_2)

object.list <- list(`NKG2A+NK` = cellchat_CT, `NKG2C+NK` = cellchat_CS)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), 
                          cell.prefix = T)

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


#### Title: Customization (CellChat) ####
#' Author: Han Zhifa
#' Date: 2021.10.08
#' 

## Part I: Differential number of interactions or interaction strength among different cell types ####
print( levels (cellchat@idents$joint) )
# 
cell1 <- "CD14+Mono"; cell2 <- ""; cell3 <- "NK"
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
celltype <- "NK" # choose one celltype for visualization
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, 
                                            idents.use = celltype)#, 
                                            #signaling.exclude = "IGTB2")
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
pos.dataset = "NKG2C+NK"
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
net.up <- subsetCommunication(cellchat, net = net, datasets = "NKG2C+NK",
                              ligand.logFC = 0.2, 
                              receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, 
                                datasets = "NKG2C+NK",
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
                                levels = c("CS", "CT")) # set factor level
plotGeneExpression(cellchat, signaling = pathways.show, 
                   split.by = "datasets",
                   colors.ggplot = T)



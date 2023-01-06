#' Title: cell cell communication (CellChat)
#' Author: Han Zhifa
#' Date: 2021.09.25
#' 
#*####**************************************####*#
#       packages and workDir and data            #
######======================================######
.libPaths("/home/hanzf/.conda/envs/scs10x/lib/R/library")
library(CellChat)
library(ggalluvial)
try(library(Seurat, lib.loc = "/home/hanzf/.conda/envs/scs10x/lib/R/library"), 
    silent = T)
require(Seurat, lib.loc = "/opt/R/4.1.0/lib/R/library")
library(NMF)

setwd("/ssd/hanzf/research/scs10x/")

combined_all <- readRDS(commandArgs(T)[2])
Idents(combined_all) <- "celltypes"

species <- commandArgs(T)[3]

results_path <- commandArgs(T)[1]
CellChat_path_0 <- paste0(results_path, "CellChat/")
dir.create(CellChat_path_0)

if("experiment" %in% colnames(combined_all@meta.data)){
  CellChat_path_1 <- paste0(results_path, "CellChat/experiment/")
  dir.create(CellChat_path_1)
}

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

for(e in unique(combined_all$experiment)){
  CellChat_path <- paste0(paste0(CellChat_path_1, e, "/"))
  dir.create(CellChat_path, recursive = T)
  
  ##### creat cellchat object
  combined <- combined_all[, combined_all$experiment == e ]
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
  
  ##### cell-cell communication analysis
  for(i in c("ALL")){ #, unique(CellChatDB$interaction$annotation)
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
    cat("input the number of patterns of outgoing:")
    nPatterns <- as.numeric(readLines("stdin", n = 1))
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
    cat("input the number of patterns of incoming:")
    nPatterns <- as.numeric(readLines("stdin", n = 1))
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
  }
}


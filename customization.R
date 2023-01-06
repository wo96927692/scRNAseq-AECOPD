#### packages and directory ####
BiocManager::install("robustbase")
library(Seurat)
require(biomaRt)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)
library(ggforce)
setwd("/ssd/hanzf/research/scs10x/")

#### data ####
combined <- readRDS("/data/hanzf/run/scs10x/results_rat/cluster/combined_celltypes.rds")
combined <- readRDS("results/cluster/combined_celltypes.rds")
sce_qc <- readRDS("resultsnodoubletrank/QC/sce_qc.rds")
combined_qc <- as.Seurat(sce_qc)

#### UMAP & TSNE ####
new.cluster.ids <- c("Endo-1", 
                     "pDC", "Neutrophil-1", "B",              
                     "NK", "CD4 T", "cDC 1", "cDC 2", 
                     "CD8 T", "Neutrophil-2", "AT2", "Alveolar Mφ-1", 
                     "Neutrophil-1", "Alveolar Mφ-2", "NK", "Adventitial FB", 
                     "Endo-2", "Endo-3", "SM", "Pericyte", 
                     "Alveolar FB", "interstitial Mφ", "CD4 T", "Neuroendocrine", 
                     "proliferating NK", "B_plasmocyte-1", "cDC 2", "Mast", 
                     "Ciliated", "AT1", "B_plasmocyte-2", "Goblet", 
                     "cDC-4", "proliferating Mφ")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
# color
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

DimPlot(combined, reduction = "umap", #split.by = "orig.ident", 
        label = T, pt.size = 0.5) +
  scale_color_manual(values = allcolour) + 
  geom_mark_hull(aes(label = combined$seurat_clusters, 
                     x = combined@reductions$umap@cell.embeddings[, "UMAP_1"], 
                     y = combined@reductions$umap@cell.embeddings[, "UMAP_2"], 
                     fill = combined$seurat_clusters), 
                 expand = unit(3, "mm"), concavity = 722) +
  theme_void()
DimPlot(combined, reduction = "tsne", 
        label = T, pt.size = 0.5)

combined <- combined[, combined$seurat_clusters != 19 & 
                       combined$seurat_clusters != 23 & 
                       combined$seurat_clusters != 25 & 
                       combined$seurat_clusters != 30 & 
                       combined$seurat_clusters != 32]

table(combined$celltypes)

#### Heatmap ####


#### Gene plot ####
DefaultAssay(combined) <- "RNA"

genes_query <- read.csv("results/customization/markers of currently known NK cell types.csv")
genes_query <- rownames(combined)[grep( "IFN", rownames(combined) )]
genes_query <- c("TMEM173", "IRF3", "JAK1", "TBK1", "IFNA5", "IFNG")

genes_query %in% rownames(combined_qc)
genes_query %in% rownames(combined)
DotPlot(combined, features = genes_query[3], group.by = "celltypes", 
        scale = F, cols = c("lightblue","red")) + 
  RotatedAxis()
DotPlot(combined[, combined$celltypes %in% c("NK", "NKG2C+NK")], 
        features = genes_query, group.by = "celltypes", 
        cols = c("lightblue","red")) + 
  RotatedAxis()
DotPlot(combined[, combined$celltypes %in% c("neutrophil-1", "neutrophil-2")], 
        features = unique(genes_query), split.by = "celltypes", 
        group.by = "experiment", cols = 1:9) + 
  RotatedAxis()
DotPlot(combined[, combined$seurat_clusters %in% c("3")], 
        features = unique(genes_query), group.by = "experiment") + 
  RotatedAxis()

combined$level_define <- paste(combined$celltypes, combined$experiment, 
                               sep = "_")
DotPlot(combined, features = unique(genes_query), group.by = "level_define") + 
  RotatedAxis() + scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0) + theme_classic()+ 
  theme(axis.text.x = element_text(angle = -45,hjust = -0.1,vjust = 1)) 


combined$experiment <- factor(combined$experiment, levels = 
                                c("NC_6m", "COPD_6m"))
VlnPlot(combined[, (1:ncol(combined))], features = genes_query[4], 
        group.by = "celltypes",split.by = "experiment",  pt.size = 0)
VlnPlot(combined[, (1:ncol(combined))[combined$celltypes %in% c("MONO_2")]], 
        features = genes_query, group.by = "experiment", pt.size = 0)
VlnPlot(combined[, combined$celltypes %in% c("NK", "NKG2C+NK")], 
        features = "CD56", #group.by = "celltypes", 
        split.by = "experiment", pt.size = 0)

FeaturePlot(combined, features = genes_query, cells = (1:ncol(combined)), 
            split.by = "experiment")
FeaturePlot(combined[, #combined$celltypes %in% c("neutrophil-1", "neutrophil-2") & 
                       #combined$experiment %in% c("NC", "COPD") & 
                       combined$seurat_clusters %in% c("47", "8", "15", "43")
                     ], features = genes_query, split.by = "experiment", 
            cells = (1:ncol(combined[, #combined$celltypes %in% c("neutrophil-1", "neutrophil-2") & 
                                     #combined$experiment %in% c("NC", "COPD") & 
                                     combined$seurat_clusters %in% c("47")
                                     ])) )

RidgePlot(combined[, (1:ncol(combined))[combined$seurat_clusters %in% c("3")]], 
          features = genes_query)

FindMarkers(combined, ident.1)

## heatmap
library(dplyr)
library(patchwork)
all.markers <- read.csv("results/cluster/markers/all_markers.csv")
all.markers <- all.markers[all.markers$cluster != 19 & 
                             all.markers$cluster != 23 & 
                             all.markers$cluster != 25 & 
                             all.markers$cluster != 30 & 
                             all.markers$cluster != 32, ]
top10 <- all.markers %>% group_by(cluster) %>% 
  top_n(n = 5, wt = avg_logFC)
DoHeatmap(combined, features=top10$gene) + NoLegend()

#### barplot for celltype percentage ####
library(reshape2)
#
table(combined$experiment)
#
ggplot(combined@meta.data[, c("celltypes", "experiment")], 
       aes(x = factor(experiment, levels = rev(unique(experiment))) )) + 
  geom_bar(aes(fill = celltypes, 
               color = 1), width = .55, position = "fill") + 
  labs(x = "experiment") + NoLegend() +
  RotatedAxis() + theme_bw() + coord_flip() + 
  scale_fill_manual(values = allcolour)
#
data <- combined@meta.data[, c("celltypes", "experiment")]
data <- as.data.frame(as.matrix(table(data)))
data <- dcast(data, formula = celltypes ~ experiment)
data$AECOPD <- round(data$/sum(data$AECOPD) * 10000)
data$`Stable COPD` <- round(data$`Stable COPD`/sum(data$`Stable COPD`) * 10000)
pdata <- data.frame(celltypes = character(), 
                    experiment = character())
for(i in data$celltypes){
  for(j in colnames(data)[-1]){
    r <- data$celltypes == i
    c <- colnames(data) == j
    pdata <- rbind(pdata, data.frame(celltypes = rep(i, data[which(r), which(c)]), 
                                     experiment = rep(j, data[which(r), which(c)])))
  }
}
pdata$experiment <- factor(pdata$experiment, 
                           levels = rev(c("AECOPD", "Stable COPD")))
ggplot(pdata, 
       aes(x = factor(celltypes, levels = unique(celltypes)))) + 
  geom_bar(aes(fill = experiment, 
               color = 1), width = .75, position = "fill") + 
  labs(x = "celltypes") + coord_flip() + theme_bw()

## barplot
length(pdata$celltypes)
# plot
ggbarplot(data.frame(t(table(pdata)))[
  data.frame(t(table(pdata)))$celltypes == "NK 2", ], 
          x = "experiment", y = "Freq", 
  ylab = "NK 2", fill = "steelblue", 
          add = c("mean_se")) #+ 
  RotatedAxis() #+ # Add pairwise comparisons p-value
#  stat_compare_means(label.y = 50)                  # Add global p-value

## gene expression
table(combined$seurat_clusters)

DimPlot(combined, reduction = "umap", 
        group.by = "batch")

#### DE ####
## Differentinial expression among celltypes/clusters and enrichment analysis
# results path
#result_path <- "/data2/hanzf/run/scs10x/results_mus/results/cluster/interested/"
result_path <- "results/cluster/interested/"
dir.create(result_path)
# clusters
Idents(combined) <- combined$seurat_clusters
DEAC <- FindMarkers(combined, ident.1 = c(3), ident.2 = c(7), 
                    min.pct = 0.25, logfc.threshold = 0.25)
# celltypes
Idents(combined) <- combined$celltypes
DEAC <- FindMarkers(combined, ident.1 = c("ActMono"), ident.2 = c("CD14+Mono"), 
                    min.pct = 0.00, logfc.threshold = 0.00, 
                    assay = "regulonAUC", slot = "data", )
# experiments
Idents(combined) <- combined$experiment
DEAC <- FindMarkers(combined[, combined$celltypes %in% c("c4")], 
                    ident.1 = c("COPD"), ident.2 = c("NC"), 
                    min.pct = 0.25, logfc.threshold = 0.25)

# ****
subcombined <- combined[, combined$celltypes %in% c("NK", "proliferating NK")]
Idents(subcombined) <- subcombined$experiment
DEAC <- FindMarkers(subcombined, ident.1 = c("COPD"), 
                    ident.2 = c("NC"), 
                    min.pct = 0.25, logfc.threshold = 0.25)
# write
write.csv(DEAC, paste0(result_path, "DEG.csv"))
# enrichment analysis
genes <- rownames(DEAC)
#
species <- "human" # human | rat | mouse
#
type <- "SYMBOL"
#
OrgDb <- paste("org", switch(species, rat = "Rn", human = "Hs", mouse = "Mm"), 
               "eg", "db", sep = ".")
#
if(OrgDb %in% library()$results[,1]){}else{BiocManager::install(OrgDb)}
require(OrgDb, character.only = T)
#
genes <- bitr(genes, fromType = type, toType = "ENTREZID", OrgDb = OrgDb)
#
organism <- switch(species, rat = "rno", human = "hsa", mouse = "mmu")
# GO
for(j in c("ALL", "BP", "CC", "MF")){
  #
  ont <- j # "ALL" | "BP" | "CC" | "MF"
  
  ## GO
  ego <- enrichGO(genes$ENTREZID, 
                  OrgDb = OrgDb, 
                  ont = ont, readable = T, 
                  pvalueCutoff = 1, qvalueCutoff = 1)
  if(!is.null(ego)){
    pdf(paste0(result_path, "/EnrichmentGO_", ont, "_dot.pdf"), 
        width = 15)
    print(dotplot(ego, showCategory = 30, title="EnrichmentGO"))
    dev.off()
    
    pdf(paste0(result_path, "/EnrichmentGO_", ont, "_bar.pdf"), 
        width = 15)
    print(barplot(ego, showCategory=30, title="EnrichmentGO"))
    dev.off()
    
    if(ont != "ALL"){
      pdf(paste0(result_path, "/EnrichmentGO_", ont, "_graph.pdf"))
      print(plotGOgraph(ego))
      dev.off()
    }
  }
}
# KEGG
ekk <- enrichKEGG(genes$ENTREZID,
                  organism = organism,
                  pvalueCutoff = 1)
if(!is.null(ekk)){
  pdf(paste0(result_path, "/EnrichmentKEGG_dot.pdf"), 
      width = 15)
  print(dotplot(ekk, title="Enrichment KEGG", showCategory = 30))
  dev.off()
  
  pdf(paste0(result_path, "/EnrichmentKEGG_bar.pdf"), 
      width = 15)
  print(barplot(ekk, title="Enrichment KEGG", showCategory = 30))
  dev.off()
}

#### 火山图 ####
#加载包
library(ggplot2)
library(ggrepel)
#读取数据
Dat <- DEAC
Dat<-read.csv('results/cluster/markers/all_markers.csv',
                header = T,stringsAsFactors = F)
#确定是上调还是下调，用于给图中点上色）
threshold_P <- 1e-50; threshold_log2FC <- 0.01
Dat$threshold = factor(ifelse(Dat$p_val_adj < threshold_P & 
                                abs(Dat$avg_log2FC) >= threshold_log2FC, 
                              ifelse(Dat$avg_log2FC>= threshold_log2FC,'Up','Down'),
                              'NoSignifi'),
                       levels=c('Up','Down','NoSignifi'))

ggplot(Dat,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+ #确定点的颜色
  geom_text_repel(
    data = Dat[Dat$p_val_adj<threshold_P & abs(Dat$avg_log2FC)> threshold_log2FC,],
    aes(label = gene),
    size = 3,
    segment.color = "black", show.legend = FALSE )+ #添加关注的点的基因名
  theme_bw()+#修改图片背景
  theme(
    legend.title = element_blank() #不显示图例标题
  )+
  ylab('-log10 (p-adj)')+ #修改y轴名称
  xlab('log2 (FoldChange)')+ #修改x轴名称
  geom_vline(xintercept=c(-threshold_log2FC,threshold_log2FC),lty=3,col="black",lwd=0.5) + #添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(threshold_P),lty=3,col="black",lwd=0.5)#添加竖线p_val_adj<0.05

#### bar plot of numbers of DEGs
Dat <- subset(Dat, cluster %in% c("1", "5", "7"))
threshold_P <- 0.05; threshold_log2FC <- 0.25
Dat$threshold = factor(ifelse(Dat$p_val_adj < threshold_P & 
                                abs(Dat$avg_log2FC) >= threshold_log2FC, 
                              ifelse(Dat$avg_log2FC>= threshold_log2FC,'Up','Down'),
                              'NoSignifi'),
                       levels=c('Up','Down','NoSignifi'))
Dat <- subset(Dat, threshold %in% c('Up','Down'))
Dat$avg_log2FC <- abs(Dat$avg_log2FC)
# Basic bar plots of means +/- se with jittered points
# compare
Dat$my_levels <- as.factor(paste(Dat$cluster, Dat$threshold, sep = "_"))
Dat$my_levels <- factor(Dat$my_levels, levels = levels(Dat$my_levels)[c(6,2,4,5,1,3)])
my_comparisons <- list( levels(Dat$my_levels)[c(1,2)], levels(Dat$my_levels)[c(2,3)], 
                        levels(Dat$my_levels)[c(1,3)], levels(Dat$my_levels)[c(4,5)], 
                        levels(Dat$my_levels)[c(5,6)], levels(Dat$my_levels)[c(4,6)])
#colors <- c("red", "green", "blue", "grey")
colors <- c(rgb(255, 255, 255, maxColorValue = 256), # white
            rgb(76, 151, 210, maxColorValue = 256), # blue
            rgb(213, 215, 212, maxColorValue = 256), # gley
            rgb(238, 57, 54, maxColorValue = 256) # red
)
#my_color <- colors[1:length(unique(Ctk_PL$experiment)) + 2]
my_color <- colors[c(4,3,2,4,3,2)]

# plot 
pdf(paste0("results/customization/", "allCD14Mono_", "number of genes.pdf"), height = 4.5, 
    width = (0.7+0.5*length(unique(Dat$my_levels))) )
# Concentrate
for (i in paste(colnames(Dat)[3]) ){
  print(
    ggbarplot(Dat, x = "my_levels", y = i, 
              fill = my_color, add = c("mean_se", "jitter"), 
              title = paste("N(UP ActMono)=", nrow(Dat[Dat$my_levels == "7_Up",]),
                            "\tN(UP CD14+Mono)=", nrow(Dat[Dat$my_levels == "1_Up",]),
                            "\nN(UP TolMono)=", nrow(Dat[Dat$my_levels == "5_Up",]),
                            "\tN(DOWN ActMono)=", nrow(Dat[Dat$my_levels == "7_Down",]),
                            "\nN(DOWN CD14+Mono)=", nrow(Dat[Dat$my_levels == "1_Down",]),
                            "\tN(DOWN TolMono)=", nrow(Dat[Dat$my_levels == "5_Down",]),
                            sep = " "), 
              ylab = paste0("|", i, "|")) + guides(shape = T) + 
      stat_compare_means(comparisons = my_comparisons, label = "p.signif", #+ # Add pairwise comparisons p-value
                         method = "wilcox.test") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    #  stat_compare_means(label.y = 50)                  # Add global p-value
  )
}
dev.off()

#####

#*####**************************************####*#
#                  CellChat                      #
######======================================######
library(CellChat)
library(ggalluvial)

cellchat <- readRDS("results/cluster/CellChat/cellchat_Secreted Signaling.rds")
CellChat_path <- "results/cluster/CellChat/"
dir.create(path)

##### creat cellchat object
data.input <- combined@assays$RNA@data
identity = data.frame(group = combined$celltypes, 
                      row.names = names(combined$celltypes))
print(unique(identity$group)) # check the cell labels
cellchat <- createCellChat(data = data.input)

cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse
# Show the structure of the database
#showDatabaseCategory(CellChatDB)
#dplyr::glimpse(CellChatDB$interaction)

for(i in unique(CellChatDB$interaction$annotation)){
  CellChatDB.use <- subsetDB(CellChatDB, search = i)# "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
  cellchat@DB <- CellChatDB.use # set the used database in the object
  
  subcellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
  #future::plan("multiprocess", workers = 1) # do parallel
  subcellchat <- identifyOverExpressedGenes(subcellchat)
  subcellchat <- identifyOverExpressedInteractions(subcellchat)
  subcellchat <- projectData(subcellchat, PPI.mouse)
  ## Inference of cell-cell communication network
  subcellchat <- computeCommunProb(subcellchat)
  subcellchat <- computeCommunProbPathway(subcellchat)
  subcellchat <- aggregateNet(subcellchat)
  ## visualization
  #pathways.show <- c("TGFb") 
  #vertex.receiver = seq(1,9) # a numeric vector
  # Hierarchy plot
  #netVisual_aggregate(subcellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)
  # Circle plot
  #netVisual_aggregate(subcellchat, signaling = pathways.show, layout = "circle", vertex.size = groupSize)
  # 2
  #netAnalysis_contribution(subcellchat, signaling = pathways.show)
  subcellchat <- netAnalysis_signalingRole(subcellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  #netVisual_signalingRole(subcellchat, signaling = pathways.show)
  
  ##  
  nPatterns = 5
  pdf(file = paste0(CellChat_path, "heatmap_", i, "_incoming.pdf"))
  subcellchat <- identifyCommunicationPatterns(subcellchat, pattern = "outgoing", 
                                               k = nPatterns)
  dev.off()
  # river plot
  pdf(file = paste0(CellChat_path, "river_", i, "_outgoing.pdf"))
  netAnalysis_river(subcellchat, pattern = "outgoing")
  dev.off()
  # dot plot
  pdf(file = paste0(CellChat_path, "dot_", i, "_outgoing.pdf"))
  netAnalysis_dot(subcellchat, pattern = "outgoing")
  dev.off()
  
  ## 
  pdf(file = paste0(CellChat_path, "heatmap_", i, "_incoming.pdf"))
  subcellchat <- identifyCommunicationPatterns(subcellchat, pattern = "incoming", 
                                               k = nPatterns)
  dev.off()
  # river plot
  pdf(file = paste0(CellChat_path, "river_", i, "_incoming.pdf"))
  netAnalysis_river(subcellchat, pattern = "incoming")
  dev.off()
  # dot plot
  pdf(file = paste0(CellChat_path, "dot_", i, "_incoming.pdf"))
  netAnalysis_dot(subcellchat, pattern = "incoming")
  dev.off()
  
  subcellchat <- computeNetSimilarity(subcellchat, type = "functional")
  subcellchat <- netEmbedding(subcellchat, type = "functional")
  subcellchat <- netClustering(subcellchat, type = "functional")
  # Visualization in 2D-space
  pdf(file = paste0(CellChat_path, "NetSimilarity_functional.pdf"))
  print(netVisual_embedding(subcellchat, type = "functional"))
  print(netVisual_embeddingZoomIn(subcellchat, type = "functional"))
  dev.off()
  
  subcellchat <- computeNetSimilarity(subcellchat, type = "structural")
  subcellchat <- netEmbedding(subcellchat, type = "structural")
  subcellchat <- netClustering(subcellchat, type = "structural")
  # Visualization in 2D-space
  pdf(file = paste0(CellChat_path, "NetSimilarity_structural.pdf"))
  netVisual_embedding(subcellchat, type = "structural")
  netVisual_embeddingZoomIn(subcellchat, type = "structural")
  dev.off()
  
  ## save data
  saveRDS(subcellchat, file = paste0(CellChat_path, "cellchat_", i, ".rds"))
}

subcellchat <- readRDS(paste0(CellChat_path, 
                              "cellchat_Secreted Signaling.rds"))

groupSize <- as.numeric(table(cellchat@idents))
pathways.show <- c("CCL")
vertex.receiver = seq(1,9) # a numeric vector
# Hierarchy plot
netVisual_aggregate(subcellchat, signaling = pathways.show, vertex.receiver = vertex.receiver, vertex.size = groupSize)
# Circle plot
netVisual_aggregate(subcellchat, signaling = pathways.show, layout = "circle", vertex.size = groupSize)
# 2
netAnalysis_contribution(subcellchat, signaling = pathways.show)
subcellchat <- netAnalysis_signalingRole(subcellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netVisual_signalingRole(subcellchat, signaling = pathways.show)

##  
#nPatterns = 5
#pdf(file = paste0(CellChat_path, "heatmap_", i, "_incoming.pdf"))
#cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", 
#                                             k = nPatterns)
#dev.off()
# river plot
#pdf(file = paste0(CellChat_path, "river_", i, "_outgoing.pdf"))
#netAnalysis_river(cellchat, pattern = "outgoing")
#dev.off()
# dot plot
#pdf(file = paste0(CellChat_path, "dot_", i, "_outgoing.pdf"))
#netAnalysis_dot(cellchat, pattern = "outgoing")
#dev.off()

## 
#pdf(file = paste0(CellChat_path, "heatmap_", i, "_incoming.pdf"))
#cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", 
#                                             k = nPatterns)
#dev.off()
# river plot
#pdf(file = paste0(CellChat_path, "river_", i, "_incoming.pdf"))
#netAnalysis_river(cellchat, pattern = "incoming")
#dev.off()
# dot plot
#pdf(file = paste0(CellChat_path, "dot_", i, "_incoming.pdf"))
#netAnalysis_dot(cellchat, pattern = "incoming")
#dev.off()

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
# Visualization in 2D-space
netVisual_embedding(subcellchat, type = "functional")
netVisual_embeddingZoomIn(subcellchat, type = "functional")

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural")
netVisual_embeddingZoomIn(cellchat, type = "structural")



#'Title: SCENIC visualization
#'Author: Han Zhifa
#'Date: 2021/12/12
#'
#### package and dir ####
## packages
# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(Seurat)

# For some of the plots:
#library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)

## dir
setwd("/ssd/hanzf/research/scs10x/")
vsnDir <- "results/cluster/SCENIC/"

#### data ####
scenicLoomPath <- file.path(paste0(vsnDir, "combined_SCENIC.loom"))
motifEnrichmentFile <- file.path(paste0(vsnDir, "reg.csv"))
file.exists(scenicLoomPath)
file.exists(motifEnrichmentFile)

## Loading results from a .loom file
loom <- open_loom(scenicLoomPath)
# Read information from loom file:
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonAucThresholds_tem <- get_regulon_thresholds(loom)
regulonAucThresholds <- names(regulonAucThresholds_tem)
names(regulonAucThresholds) <- regulonAucThresholds_tem
embeddings <- get_embeddings(loom)
#cellClusters <- get_clusterings(loom)
#cellClusters <- subset(cellClusters, select = "celltypes")
if (T) {
  combined <- readRDS("results/cluster/combined_celltypes.rds")
  cellClusters <- data.frame(celltypes = combined$celltypes)
}
close_loom(loom)
# the motif enrichment results
motifEnrichment <- data.table::fread(motifEnrichmentFile, header=T, skip=1)[-3,]
colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")

#### Scoring the network activity ####
selectedResolution <- "celltypes" # select resolution
# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellClusters), cellClusters[,selectedResolution]) 
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
# plot:
#options(repr.plot.width=8, repr.plot.height=10) # To set the figure size in Jupyter
pdf(paste0(vsnDir, "Heatmap_regulonActivity.pdf"))
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[, -11], name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=1))) # row font size
dev.off()
regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators)
viewTable(topRegulators, options = list(pageLength = 10))
write.csv(topRegulators, file = paste0(vsnDir, "topRegulators.csv"))

#### Regulon Specificity Score (RSS) ####
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellClusters[colnames(regulonAUC), 
                                                                   selectedResolution])
write.csv(topRegulators, file = paste0(vsnDir, "RSS.csv"))
# Showing regulons and cell types with any RSS > 0.01 
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

#options(repr.plot.width=2, repr.plot.height=2) # To set the figure size in Jupyter
#plotRSS_oneSet(rss, setName = "NEU") # cluster ID (celltype name)
pdf(paste0(vsnDir, "celltype specific regulons.pdf"), 
    height = 4, width = 4)
for ( i in unique(as.character(cellClusters$celltypes)) ) { print(i)
  print( plotRSS_oneSet(rss, setName = i) ) # cluster ID (celltype name)
}
dev.off()

#### Cell types based on the GRN activity ####
cat(names(embeddings), sep="\n")
# Overview of these embeddings (see below for details)
regulonsToPlot <- "EGR2"
#options(repr.plot.width=10, repr.plot.height=8) # To set the figure size in Jupyter
par(mfrow=c(2, ceiling(length(names(embeddings))/2)))
for (selectedEmbedding in names(embeddings)){
  AUCell::AUCell_plotTSNE(embeddings[[selectedEmbedding]], exprMat_log, 
                          regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, 
                          sub=selectedEmbedding)
}
  
selectedEmbedding <- embeddings[["umap.rna"]] # change if desired...
tfsToPlot <- c("ETS1", "SPI1", "EBF1") 
regulonsToPlot <- unlist(lapply(tfsToPlot, 
                                function(x) grep(paste0("^", x), 
                                                 rownames(regulonAUC), value=TRUE)))

#options(repr.plot.width=10, repr.plot.height=8) # To set the figure size in Jupyter
par(mfrow=c(2,3))
# Plot expression:
AUCell::AUCell_plotTSNE(selectedEmbedding, exprMat_log[tfsToPlot,], plots=c("Expression"), cex = .5)
# Plot regulon activity:
AUCell::AUCell_plotTSNE(selectedEmbedding, exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5)

# Subset some cells to make the plots faster:
nCells <- 1000
set.seed(123)
cellsSelected <- sample(colnames(regulonAUC), nCells) 
regulonAUC_subset <- regulonAUC[regulonOrder, which(colnames(regulonAUC) %in% cellsSelected)]
dim(regulonAUC_subset)
selectedEmbedding_subset <- selectedEmbedding[colnames(regulonAUC_subset), ]

# Save AUC as PDF:
pdf(paste0(vsnDir, "RegulonActivity.pdf"), width=20, height=15)
par(mfrow=c(4,6))
AUCell::AUCell_plotTSNE(selectedEmbedding_subset, 
                        cellsAUC=regulonAUC_subset, 
                        plots="AUC")
dev.off()

#options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(selectedEmbedding, .5)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=6, drawlabels=FALSE)

regulonsToPlot <- "EGR2(+)"
#options(repr.plot.width=10, repr.plot.height=4) # To set the figure size in Jupyter
par(mfrow=c(1,3))
AUCell::AUCell_plotTSNE(selectedEmbedding, exprMat_log, 
                        regulonAUC[regulonsToPlot,], 
                        thresholds = regulonAucThresholds[regulonsToPlot],
                        plots=c("AUC", "histogram", "binary"), cex = .5)
# shiny
aucellApp <- AUCell_createViewerApp(auc=regulonAUC,
                                    thresholds=regulonAucThresholds,
                                    tSNE=selectedEmbedding, 
                                    exprMat=exprMat_log)
savedSelections <- shiny::runApp(aucellApp)

#### binarizeAUC ####
# This function will be included in the next version of AUCell
if (T) {
  binarizeAUC <- function(auc, thresholds)
  {
    thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
    regulonsCells <- setNames(lapply(names(thresholds), 
                                     function(x) {
                                       trh <- thresholds[x]
                                       names(which(getAUC(auc)[x,]>trh))
                                     }),names(thresholds))
    
    regulonActivity <- reshape2::melt(regulonsCells)
    binaryRegulonActivity <- 
      t(table(regulonActivity[,1], regulonActivity[,2]))
    class(binaryRegulonActivity) <- "matrix"
    
    return(binaryRegulonActivity)
  }
}

binaryRegulonActivity <- 
  binarizeAUC(regulonAUC, regulonAucThresholds)
dim(binaryRegulonActivity)
binaryRegulonActivity[1:5,1:3]

# Subset some cells to make the plots faster:
nCells <- 1000
set.seed(123)
cellsSelected <- sample(colnames(regulonAUC), nCells)
binAct_subset <- 
  binaryRegulonActivity[, which(colnames(binaryRegulonActivity) %in% cellsSelected)]
dim(binAct_subset)

options(repr.plot.width=12, repr.plot.height=10) # To set the figure size in Jupyter
binAct_subset <- binAct_subset[regulonOrder,]
sub_cellClusters <- subset(cellClusters, 
                           rownames(cellClusters) %in% colnames(binAct_subset))
my_colors <- list(celltypes = c(`CD8+T` = "#4DA0DB", `CD14+Mono` = "#DB5A58", 
                                `NK` = "#53B7DB", `CD4+T` = "#DB5A37",
                                `NEU(1)` = "#37DB6B", `TolMono` = "#DB942C",
                                `B(1)` = "#DBCE3D", ActMono = "#DB3B9C",
                                `B(2)` = "#DB8A32",`NKG2C+NK` = "#DB5B97",
                                `CD16+Mono` = "#459ADB", `DC` = "#42DBC9",
                                `NEU(2)` = "#DB9E42", Platelet = "#078F67",
                                Plasma = "#15708F"))
  
ComplexHeatmap::pheatmap(binAct_subset, name="Binarized activity", 
                        col = c("white", "black"), 
                        annotation = T,
                        annotation_col = sub_cellClusters, 
                        #annotation_colors = my_colors,
                        legend = F, show_rownames = T, 
                        fontsize_row = 1.5, 
                        show_colnames = F
                        ) # row font size

#### combine Seurat and loom ####
combined <- readRDS("results/cluster/combined_celltypes.rds")
combined[["regulonAUC"]] <- CreateAssayObject( 
  data = as.matrix(regulonAUC@assays@data$AUC) )
rss <- data.frame(rss)
rss$regulor <- rownames(rss)
# create dir
DEAUC_path <- paste0(paste0(vsnDir, "DEregulons/"))
dir.create(DEAUC_path, recursive = T)
for(i in unique (combined$celltypes) ){
  # DE
  DEs <- FindMarkers(combined, ident.1 = "AECOPD", ident.2 = "Stable COPD", 
                     group.by = "experiment", subset.ident = i, 
                     assay = "regulonAUC", slot = "data", 
                     logfc.threshold = 0)
  # output
  write.csv(DEs, 
            file= paste0(DEAUC_path, i, ".csv"), 
            quote=F, row.names=T)
  
  # plot
  Dat <- DEs; Dat$regulor <- rownames(Dat)
  if( all(rownames(Dat) %in% rownames(rss)) ){
    i2 = data.frame(a = 1); colnames(i2) = i
    Dat <- merge(Dat, rss[c(gsub("[\\+\\(\\)]", "\\.", i), "regulor")], 
                 by = "regulor", all.y = F)
    rownames(Dat) <- Dat$regulor
  }
  Dat$threshold = factor(ifelse(Dat$p_val_adj < 0.05 & 
                                  abs(Dat$avg_log2FC) >= 0.02, 
                                ifelse(Dat$avg_log2FC>= 0.02,'Up','Down'),
                                'NoSignifi'),
                         levels=c('Up','Down','NoSignifi'))
  pdf(file = paste0(DEAUC_path, i, ".pdf"), height = 5, width = 6)
    print(
    ggplot( Dat,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold) )+
      geom_point(size = 100^((Dat[[gsub("[\\+\\(\\)]", "\\.", i)]])) )+
      scale_color_manual(values=c("#DC143C","#00008B","#808080"))+ #确定点的颜色
      geom_text_repel(
        data = Dat[Dat$p_val_adj<0.05 & abs(Dat$avg_log2FC)> 0.02,],
        aes(label = regulor),
        size = 3,
        segment.color = "black", show.legend = FALSE )+ #添加关注的点的基因名
      theme_bw()+#修改图片背景
      theme(
        legend.title = element_blank() #不显示图例标题
      )+
      ylab('-log10 (p-adj)')+ #修改y轴名称
      xlab('log2 (FoldChange of regulor AUC)')+ #修改x轴名称
      geom_vline(xintercept=c(-0.02,0.02),lty=3,col="black",lwd=0.5) + #添加横线|FoldChange|>2
      geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)#添加竖线p_val_adj<0.05
    )
  dev.off()
}


#'10x data QC
#'Data: 2021/09/25
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
library(SingleCellExperiment)
require(scater)
try(library(Seurat, lib.loc = "/opt/R/4.1.0/lib/R/library"), 
    silent = T)
require(Seurat, lib.loc = "/home/hanzf/.conda/envs/scs10x/lib/R/library")


#getwd()
setwd("/ssd/hanzf/research/scs10x")

#####

#*####**************************************####*#
#                    get data                    #
######======================================######
## read data
sample_name <- commandArgs(T)[1]
sample_name <- unlist(strsplit(sample_name, ","))
data_path <- commandArgs(T)[2]
cluster_path <- commandArgs(T)[3]
dir.create(cluster_path)
results_path <- "results/QC/"
dir.create(results_path)
# The first col of annotation file must be sample name as same as sample_name
annotation <- read.csv(paste0(data_path, "sample_annotation.csv"), 
                       row.names = 1, stringsAsFactors = F)
if(identical(rownames(annotation), sample_name)){}else{
  cat("Please note the sample annotation file have some problem")
}
# add annotation and doublet information
counts <- list()
for(i in 1:length(sample_name)){
  counts[[i]] <- Read10X(paste0(data_path, sample_name[i], "/outs/filtered_feature_bc_matrix/"))
  colnames(counts[[i]]) <- paste(sample_name[i], colnames(counts[[i]]), sep = "-")
  names(counts)[i] <- sample_name[i]
  counts[[i]] <- CreateSeuratObject(counts = counts[[i]], assay = "RNA", 
                                    project = sample_name[i])
  #counts[[i]]@meta.data <- cbind(counts[[i]]@meta.data, as.data.frame(sapply(annotation, 
  #                                                            rep, ncol(counts[[i]]))))
  addmetadata <- as.data.frame(sapply(annotation[i, ], rep, ncol(counts[[i]])))
  ## read doublet scores
  doublet_scores <- read.table(paste0("results/doublet/", sample_name[i], "_doublet_scores.txt"), 
                               header = F, stringsAsFactors = F)
  colnames(doublet_scores) <- "doublet_scores"
  doublet <- read.table(paste0("results/doublet/", sample_name[i], "_predicted_doublets.txt"), 
                        header = F, stringsAsFactors = F)
  colnames(doublet) <- "doublet"
  doublet$doublet <- as.logical(doublet$doublet)
  ## add the annotation information
  addmetadata <- as.data.frame(sapply(annotation[i, ], rep, ncol(counts[[i]])))
  rownames(addmetadata) <- colnames(counts[[i]])
  addmetadata <- cbind(addmetadata, doublet_scores, data.frame(doublet = doublet))
  counts[[i]] <- AddMetaData(counts[[i]], addmetadata)
}

## combine all seurat data
all_counts <- merge(counts[[1]], counts[-1])

## map seurat data to sce
# define function
#Seurat2SCE <- function(SeuratObject){
#  result <- SingleCellExperiment(assays = list(counts = SeuratObject@assays$RNA@counts), 
#                                 colData = SeuratObject@meta.data, 
#                                 rowData = SeuratObject@assays$RNA@meta.features)
#}
# transform
#sce <- Seurat2SCE(all_counts)
sce <- as.SingleCellExperiment(all_counts)

## identify mitochondria genes
MT.logi <- grepl(pattern = "^[Mm][tT]-", x = rownames(sce))
all_counts[["percent.mt"]] <- PercentageFeatureSet(all_counts, pattern = "^[Mm][tT]-")

#####

#*####**************************************####*#
#                    cell QC                     #
######======================================######
## conduct "scater" QC metrics
per.cell <- perCellQCMetrics(sce, subset = list(MT = MT.logi))
sce_qc <- sce
colData(sce_qc) <- cbind(colData(sce_qc), per.cell)

## visualization
# detected genes vs sum counts
png(paste0(results_path, "genes_counts.png"), 
    width = 480*2, height = 480*2)
plotColData(sce_qc, x = "sum", y="detected", colour_by="orig.ident")
dev.off()
# total_counts
pdf(paste0(results_path, "hist_total_counts.pdf"))
hist(sce_qc$sum, breaks = 500, xlim = c(0, 35000))
abline(v=(median(sce_qc$sum) + 3*mad(sce_qc$sum)), col = "red")
abline(v=(median(sce_qc$sum) - 3*mad(sce_qc$sum)), col = "blue")
dev.off()
# total_genes_by_counts
pdf(paste0(results_path, "hist_total_features_by_counts.pdf"))
hist(sce_qc$detected, breaks = 200)
abline(v=(median(sce_qc$detected) + 3*mad(sce_qc$detected)), col="red")
abline(v=(median(sce_qc$detected) - 3*mad(sce_qc$detected)), col="blue")
dev.off()
# MT percent
pdf(paste0(results_path, "hist_pct_counts_MT.pdf"))
hist(sce_qc$subsets_MT_percent, breaks = 50)
abline(v=(median(sce_qc$subsets_MT_percent) + 3*mad(sce_qc$subsets_MT_percent)), col="red")
abline(v=(median(sce_qc$subsets_MT_percent) - 3*mad(sce_qc$subsets_MT_percent)), col="blue")
dev.off()
png(paste0(results_path, "dot_MTvsFeatures.png"), 
    width = 480*2, height = 480*2)
plotColData(sce_qc, x = "detected", y = "subsets_MT_percent", 
            colour_by = "orig.ident", shape_by = "batch")
dev.off()
png(paste0(results_path, "dot_MTvsFeatures_facet.png"), 
    height = 350, 
    width = length(unique(sce_qc$orig.ident)) * 350)
plotColData(sce_qc, x = "sum", y="subsets_MT_percent", 
            other_fields = "orig.ident") + facet_wrap(~orig.ident)
dev.off()
# output with Seurat
png(paste0(results_path, "Seurat_qc.png"), width = 15*70)
VlnPlot(all_counts, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

## Visulizing doublets
sce_qc <- runPCA(sce_qc)

png(paste0(results_path, "PCA_doublets.png"), 
    width = 480*2, height = 480*2)
plotReducedDim(sce_qc, dimred = "PCA", 
               shape_by = "doublet", 
               colour_by = "orig.ident", 
               size_by = "doublet_scores")
dev.off()

## plot the top genes percentage of total counts
pdf(paste0(results_path, "highest_expressive_genes.pdf"))
plotHighestExprs(sce_qc, exprs_values = "counts")
dev.off()

## intergrat the several qc
cat("Please choose the appropriate thresholds according the pictures in dir: results/QC\n")
cat("max threshold of pecentage for MT: ")
threshold_pct_counts_MT <- readLines("stdin", n = 1)
print(threshold_pct_counts_MT)
cat("mix threshold of number for total features: ")
mix_threshold_total_features_by_counts <- readLines("stdin", n = 1)
print(mix_threshold_total_features_by_counts)
cat("max threshold of number for total features: ")
max_threshold_total_features_by_counts <- readLines("stdin", n = 1)
print(max_threshold_total_features_by_counts)
cat("mix threshold of number for total counts: ")
mix_threshold_total_counts <- readLines("stdin", n = 1)
print(mix_threshold_total_counts)
cat("max threshold of number for total counts: ")
max_threshold_total_counts <- readLines("stdin", n = 1)
print(max_threshold_total_counts)

colData(sce_qc)$QC <- !sce_qc$doublet & 
  sce_qc$subsets_MT_percent < as.numeric(threshold_pct_counts_MT)  & 
  sce_qc$detected > as.numeric(mix_threshold_total_features_by_counts) & 
  sce_qc$detected < as.numeric(max_threshold_total_features_by_counts) & 
  sce_qc$sum > as.numeric(mix_threshold_total_counts) & 
  sce_qc$sum < as.numeric(max_threshold_total_counts) & 
  !sce_qc$doublet
print(table(sce_qc$QC))

## quick outliner
qc.stats <- quickPerCellQC(per.cell, percent_subsets="subsets_MT_percent")

colData(sce_qc)$outlier <- qc.stats$discard

## output QC condition
qc_condition <- data.frame(threshold_pct_counts_MT, 
                           mix_threshold_total_features_by_counts, 
                           max_threshold_total_features_by_counts, 
                           mix_threshold_total_counts, 
                           max_threshold_total_counts, 
                           mix_gene_express = 2)
# the mix_gene_express used in gene QC
write.csv(qc_condition, paste0(results_path, "qc_condition.csv"), row.names = F)

## rm PCA results for call outliers
reducedDim(sce_qc) <- NULL

#####

#*####**************************************####*#
#                    gene QC                     #
######======================================######
## gene level QC
per.feat <- perFeatureQCMetrics(sce)
rowData(sce_qc) <- cbind(rowData(sce_qc), per.feat)

## obstract the genes express (counts > 0) in more than 2 cells
print("obstract the genes expressed (counts > 0) in more than 2 cells\n")
keep_feature <- nexprs(sce, byrow = T, detection_limit = 0) >= 2
table(keep_feature)

## intergrat several gene qc
rowData(sce_qc)$QC <- keep_feature
rowData(sce_qc)$MT <- MT.logi

is_endog_genes <- !(rowData(sce_qc)$MT)

#####

#*####**************************************####*#
#           Assess confounding factors           #
######======================================######
## 
# tmp is the sce_qc without endog genes
tmp <- runPCA(sce_qc[is_endog_genes, ], 
              exprs_values = "logcounts", 
              ncomponents = 100)
# PCA plot with QC as shape
pdf(paste0(results_path, "PCA_QC.pdf"))
plotPCA(tmp, colour_by = "orig.ident", size_by = "detected", 
        shape_by = "QC", ncomponents = 4)
dev.off()

## 1st PC vs variables of colData colnames
# 1st PC vs "total_features_by_counts" of colData colnames
pdf(paste0(results_path, "PCvsFeatures.pdf"))
plotExplanatoryPCs(tmp, variables = "detected")
dev.off()
# 1st PC vs "total_counts" of colData colnames
pdf(paste0(results_path, "PCvsCounts.pdf"))
plotExplanatoryPCs(tmp, variables = "sum")
dev.off()

## total correlation between each gene and variables of colData colnames
# plot 
pdf(paste0(results_path, "explanatory_variables.pdf"))
plotExplanatoryVariables(tmp, exprs_values = "logcounts", 
                         variables = c("detected", "sum", "subsets_MT_percent", 
                                       "orig.ident", colnames(annotation)))
dev.off()

#####


#*####**************************************####*#
#                   PCA/tSNE                     #
######======================================######
cat("execute PCA (tSNE) (it will take long time) (T or F): \n")
doPCA_logi <- as.logical(readLines("stdin", n = 1))
if(doPCA_logi) {
  ## check the amounts of QCed matrix
  #dim(sce_qc[rowData(sce_qc)$QC, sce_qc$QC])
  
  ## build the amounts of QCed matrix
  sce_QC <- sce_qc[rowData(sce_qc)$QC & is_endog_genes, sce_qc$QC]
  
  ## plot PCA before and after QC
  # before QC with counts
  tmp <- runPCA(sce_qc[is_endog_genes, ], exprs_values = "counts")
  reducedDimNames(tmp)
  pdf(paste0(results_path, "PCA_beforeQC.pdf"))
  print(plotPCA(tmp, shape_by = "QC", 
                colour_by = "experiment", size_by = "detected"))
  dev.off()
  # after QC with counts
  tmp <- runPCA(sce_QC, exprs_values = "counts")
  reducedDimNames(tmp)
  pdf(paste0(results_path, "PCA_afterQC.pdf"))
  print(plotPCA(tmp, colour_by = "experiment", size_by = "detected", 
                shape_by = "batch"))
  dev.off()
  
  # before QC with logcounts_raw
  tmp <- runPCA(sce_qc[is_endog_genes, ], exprs_values = "logcounts")
  reducedDimNames(tmp)
  pdf(paste0(results_path, "PCA_LOG_beforeQC.pdf"))
  print(plotPCA(tmp, shape_by = "outlier", size_by = "detected", 
                colour_by = "orig.ident"))
  dev.off()
  # after QC with counts
  tmp <- runPCA(sce_QC, exprs_values = "logcounts")
  reducedDimNames(tmp)
  pdf(paste0(results_path, "PCA_LOG_afterQC.pdf"))
  print(plotPCA(tmp, size_by = "detected", colour_by = "orig.ident", 
                shape = "batch"))
  dev.off()
  
  ## plot tSNE before and after QC based on logcounts
  # before QC with counts
  set.seed(12345)
  tmp <- runTSNE(sce_qc[is_endog_genes, ], exprs_values = "logcounts", perplexity = 50)
  reducedDimNames(tmp)
  pdf(paste0(results_path, "tSNE_LOG_beforeQC.pdf"))
  print(plotTSNE(tmp, shape_by = "QC", 
                 colour_by = "experiment", size_by = "detected"))
  dev.off()
  # after QC with counts
  tmp <- runTSNE(sce_QC, exprs_values = "logcounts", perplexity = 50)
  reducedDimNames(tmp)
  pdf(paste0(results_path, "tSNE_LOG_afterQC.pdf"))
  print(plotTSNE(tmp, colour_by = "experiment", size_by = "detected", 
                 shape_by = "batch"))
  dev.off()
}

#####


#*####**************************************####*#
#                   save data                    #
######======================================######
## save data
write.csv(colData(sce_qc), file = paste0(results_path, "/metadata_samples.csv"), row.names = T)
#saveRDS(sce, file = paste0(results_path, "/sce.rds"))
#save.image(file = paste0(results_path, "/qc.RData"))
saveRDS(sce_qc, file = paste0(results_path, "/sce_qc.rds"))
##
seurat_qc <- as.Seurat(sce_qc[rowData(sce_qc)$QC, sce_qc$QC])
saveRDS(seurat_qc, file = paste0(cluster_path, "/seurat_QC.rds"))
#save.image(paste0(results_path, "/QC.RData"))

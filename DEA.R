#'DE scs
#*####**************************************####*#
#             packages and workDir               #
######======================================######
#install.packages("jpeg")
#BiocManager::install("DESeq2")
#BiocManager::install("biomaRt")
#install.packages("pheatmap")
#install.packages("amap")
#install.packages(c("gplots", "amap", "ggplot2"))
try(library(Seurat, lib.loc = "/home/hanzf/.conda/envs/scs10x/lib/R/library"), 
    silent = T)
require(Seurat, lib.loc = "/opt/R/4.1.0/lib/R/library")
library(biomaRt)
library(ggplot2)
library(gplots)
library(amap)

## pathway analysis
#BiocManager::install(version='3.10')
#BiocManager::install("topGO")
#BiocManager::install("SingleR")
#BiocManager::install("org.Rn.eg.db")
#BiocManager::install("Rgraphviz")
#BiocManager::install("pathview")
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(dplyr, lib.loc = "/opt/R/4.1.0/lib/R/library")

# 设制路径
#getwd()
setwd("/ssd/hanzf/research/scs10x")

#*####**************************************####*#
#                    get data                    #
######======================================######
##
results_path <- paste0(commandArgs(T)[2], "DEA/")
dir.create(results_path)

## read combined
combined <- readRDS(commandArgs(T)[1])

#####

#*####**************************************####*#
#   differential experiment analysis(monocle3)   #
######======================================######
# definite cds class
DefaultAssay(combined) <- "RNA"
#Idents(combined) <- "celltypes"

## Delete the celltypes who have less than 10 cells
deleted_celltypes <- rownames(table(combined$celltypes, combined$experiment))[
  which(table(combined$celltypes, combined$experiment) < 10)
]
combined <- combined[, !(combined$celltypes %in% deleted_celltypes) ]

##
for(i in unique (combined$celltypes) ){
  ## create dir
  DEAC_path <- paste0(results_path, i, "/")
  dir.create(DEAC_path, recursive = T)
  
  ## DE
  DEG <- FindMarkers(combined, ident.1 = "AECOPD", ident.2 = "Stable COPD", 
                     group.by = "experiment", subset.ident = i)
  # output
  write.csv(DEG, 
            file= paste0(DEAC_path, "DEG.csv"), 
            quote=F, row.names=T)
  #
  sig_gene <- subset(DEG, p_val_adj < 0.05)
  
  ## enrichment analysis
  for(j in c("ALL", "BP", "CC", "MF")){
    ## get gene list and it's annotation
    
    genes <- rownames(sig_gene)
    
    species <- "human" # human | rat | mouse
    
    type <- "SYMBOL" # "SYMBOL" | "ENSEMBL"
    
    ont <- j # "ALL" | "BP" | "CC" | "MF"
    
    OrgDb <- paste("org", switch(species, rat = "Rn", human = "Hs", mouse = "Mm"), 
                   "eg", "db", sep = ".")
    
    if(OrgDb %in% library()$results[,1]){}else{BiocManager::install(OrgDb)}
    require(OrgDb, character.only = T)
    
    genes <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
    
    organism <- switch(species, rat = "rno", human = "hsa", mouse = "mmu")
    
    ## GO
    ego <- enrichGO(genes$ENTREZID, 
                    OrgDb = OrgDb, 
                    ont = ont, readable = T, 
                    pvalueCutoff = 1, qvalueCutoff = 1)
    if(!is.null(ego)){
      pdf(paste0(DEAC_path, "EnrichmentGO_", ont, "_dot.pdf"), 
          width = 15)
      print(dotplot(ego))
      dev.off()
      
      pdf(paste0(DEAC_path, "EnrichmentGO_", ont, "_bar.pdf"), 
          width = 15)
      print(barplot(ego, showCategory=20, title="EnrichmentGO"))
      dev.off()
      
      if(ont != "ALL"){
        pdf(paste0(DEAC_path, "/EnrichmentGO_", ont, "_graph.pdf"))
        print(plotGOgraph(ego))
        dev.off()
      }
    }
    
    ## KEGG
    ekk <- enrichKEGG(genes$ENTREZID,
                      organism = organism,
                      pvalueCutoff = 1)
    if(!is.null(ekk)){
      pdf(paste0(DEAC_path, "/EnrichmentKEGG_dot.pdf"), 
          width = 15)
      print(dotplot(ekk, title="Enrichment KEGG"))
      dev.off()
      
      pdf(paste0(DEAC_path, "/EnrichmentKEGG_bar.pdf"), 
          width = 15)
      print(barplot(ekk))
      dev.off()
    }
  }
}


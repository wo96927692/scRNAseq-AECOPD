#' Title: Integration and MNN
#' Author: Han Zhifa
#' Date: 2022.4.11
#' 
#*####**************************************####*#
#       packages and workDir and data            #
######======================================######
.libPaths()
#BiocManager::install(version = "3.10", lib = "/usr/lib64/R/library")
try(library(Seurat, lib.loc = "/home/hanzf/.conda/envs/scs10x/lib/R/library"), 
    silent = T)
require(Seurat, lib.loc = "/opt/R/4.1.0/lib/R/library")
results_path <- commandArgs(T)[1]
seurat <- readRDS(commandArgs(T)[2])
## remove the subjects who have less 50 cells
seurat <- seurat[, seurat$orig.ident %in% 
                   names(table(seurat$orig.ident))[table(seurat$orig.ident) > 50]]

######
cat("Do you want to integrate data?(T or F):")
wthItg <- as.logical(readLines("stdin", n = 1))
if(!wthItg){
  combined <- NormalizeData(seurat)
  ## scale data
  cat("Do you want to scale all genes?(T of F): ")
  wthAllGenes <- as.logical(readLines("stdin", n = 1))
  if(wthAllGenes){
    features = rownames(combined)
  }else{combined <- FindVariableFeatures(combined)}
  cat("The factors can be adjusted include: "); print(colnames(combined[[]]))
  cat("the usually adjusted factors include: nCount_RNA,subsets_MT_percent,batch,CC.Difference\n")
  cat("Or just input 'DonotRegress'")
  cat("do you want to define the adjusted factors(T or F): ")
  log <- as.logical(readLines("stdin", n = 1))
  if(log){
    cat("which factors you want to adjust(',' sparated):")
    factors <- unlist(strsplit(readLines("stdin", n = 1), split = ","))
    cat("\nOr just input 'DonotRegress'\n")
    if("DonotRegress" == factors) factors <- NULL
    if(wthAllGenes){
      combined <- ScaleData(combined, vars.to.regress = factors, 
                            features = features)
    }else{combined <- ScaleData(combined, vars.to.regress = factors)}
    
  }else{
    if(wthAllGenes){
      combined <- ScaleData(combined, features = features)
    }else{combined <- ScaleData(combined)}
  }
}else{
  cat("The factors can to split include (only Quality traits): ")
  print(colnames(seurat[[]]))
  cat("select one factor to split: ")
  split_factor <- readLines("stdin", n = 1)
  combined.list <- SplitObject(seurat, split.by = split_factor)
  combined.list <- lapply(X = combined.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", 
                              nfeatures = 2000)
  })
  combined.features <- SelectIntegrationFeatures(object.list = combined.list)
  #
  cat("Which reduction method you select?(cca, rpca, rlsi): ")
  reduction_method <- readLines("stdin", n = 1)
  # strength of reduction method
  rpca_strength <- 5
  if(reduction_method == "rpca"){
    cat("The default n.anchor of rpca is 5, please enter your n.anchor): ")
    rpca_strength <- as.numeric(readLines("stdin", n = 1))
  }
  
  #
  cat("Do you want to use SCTransform to integrate different list?(T or F):")
  SCTlog <- as.logical(readLines("stdin", n = 1))
  
  if(SCTlog){
    ## scale data
    cat("The factors can be adjusted include: "); print(colnames(seurat[[]]))
    cat("the default adjusted factors include: nCount_RNA,percent.mt,batch,CC.Difference")
    cat("do you want to redefine the adjusted factors(T or F): ")
    log <- as.logical(readLines("stdin", n = 1))
    if(log){
      cat("which factors you want to adjust(',' sparated):")
      factors <- unlist(strsplit(readLines("stdin", n = 1), split = ","))
      cat("\nOr just input 'DonotRegress'\n")
      if("DonotRegress" == factors) factors <- NULL
      for (i in names(combined.list)) {
        combined.list[[i]] <- SCTransform(combined.list[[i]], verbose = FALSE, 
                                          vars.to.regress = factors)
      }
    }else{
      for (i in names(combined.list)) {
        combined.list[[i]] <- SCTransform(combined.list[[i]], verbose = FALSE, 
                                          vars.to.regress = c("nCount_RNA", "percent.mt", 
                                                              "batch", "CC.Difference"))
      }
    }
    #
    combined.list <- PrepSCTIntegration(object.list = combined.list, anchor.features = combined.features)
    combined.anchors <- FindIntegrationAnchors(object.list = combined.list, normalization.method = "SCT", 
                                               anchor.features = combined.features)
    combined <- IntegrateData(anchorset = combined.anchors, normalization.method = "SCT")
    #
    combined <- RunPCA(object = combined, verbose = FALSE)
    combined <- RunUMAP(object = combined, dims = 1:30)
  }else{
    LNIntegrate <- function(counts, features, reduction_method){
      ## intergration
      if(reduction_method == "rpca"){
        counts <- lapply(X = counts, FUN = function(x) {
          x <- ScaleData(x, features = features, verbose = FALSE)
          x <- RunPCA(x, features = features, verbose = FALSE)
        })
      }
      anchors <- FindIntegrationAnchors(object.list = counts, dims = 1:30, 
                                        reduction = reduction_method, 
                                        k.anchor = rpca_strength, 
                                        anchor.features = features)
      combined <- IntegrateData(anchorset = anchors, k.weight = 30)
    }
    combined <- LNIntegrate(combined.list, combined.features, reduction_method)
    
    ## 
    DefaultAssay(combined) <- "integrated"
    
    ## scale data
    cat("The factors can be adjusted include: "); print(colnames(combined[[]]))
    cat("the default adjusted factors include: nCount_RNA,subsets_MT_percent,batch,CC.Difference")
    cat("do you want to redefine the adjusted factors(T or F): ")
    log <- as.logical(readLines("stdin", n = 1))
    if(log){
      cat("which factors you want to adjust(',' sparated):")
      factors <- unlist(strsplit(readLines("stdin", n = 1), split = ","))
      cat("\nOr just input 'DonotRegress'\n")
      if("DonotRegress" == factors) factors <- NULL
      combined <- ScaleData(combined, vars.to.regress = factors)
    }else{
      combined <- ScaleData(combined, 
                            vars.to.regress = c("nCount_RNA", "subsets_MT_percent", "batch", 
                                                "CC.Difference"))
    }
    
  }
}


#####

## save data
#file.remove(commandArgs(T)[2])
saveRDS(combined, file = paste0(results_path, "combined_integration.rds"))
if(wthItg){
  write.csv(data.frame(split_factor = split_factor), file = 
              paste0(results_path, "split_factor.csv"), row.names = F)
}else{write.csv(data.frame(split_factor = "NULL"), file = 
                  paste0(results_path, "split_factor.csv"), row.names = F)}

if("log" %in% ls() & log){
  write.csv(data.frame(regress_factors = factors), file = 
              paste0(results_path, "regress_factor.csv"), row.names = F)
}else{
  if(("SCTlog" %in% ls()) & SCTlog){ }else{
    write.csv(data.frame(regress_factors = factors), file = 
                paste0(results_path, "regress_factor.csv"), row.names = F)
  }
}

getwd()

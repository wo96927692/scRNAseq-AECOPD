#'Title: Velocity
#'Author: Han Zhifa
#'Date: 2021/1/4
#'
#### package and dir ####
## dir
setwd("/home/hanzf/Research/run/scs10x/")

## loading packages
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

library(devtools)
devtools::install_github('kharchenkolab/pagoda2', build_vignettes = TRUE)
library(pagoda2)

#### data ####
curl::curl_download(url = 'http://pklab.med.harvard.edu/velocyto/mouseBM/SCG71.loom', destfile
                    = 'logs/SCG71.loom')
                    
ldat <- ReadVelocity(file = "velocity/merged.loom")
bm <- as.Seurat(x = ldat)
bm <- SCTransform(object = bm, assay = "spliced")
bm <- RunPCA(object = bm, verbose = FALSE)
bm <- FindNeighbors(object = bm, dims = 1:20)
bm <- FindClusters(object = bm)
bm <- RunUMAP(object = bm, dims = 1:20)
bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
names(x = ident.colors) <- levels(x = bm)
cell.colors <- ident.colors[Idents(object = bm)]
names(x = ce
      ll.colors) <- colnames(x = bm)
show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), 
                               vel = Tool(object = bm, slot = "RunVelocity"), 
                               n = 200, scale = "sqrt", 
                               cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)


ldat <- read.loom.matrices("velocity/merged.loom")

emat <- ldat$spliced
emat <- emat[,colSums(emat)>=1e3]

r <- Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T)

r$adjustVariance(plot=T,do.par=T,gam.k=10)

r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)

r$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine')

r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)

par(mfrow=c(1))
r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,
                min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,
                main='cell clusters')
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"Xist"],main='Xist')


emat <- ldat$spliced; nmat <- ldat$unspliced
emat <- emat[,rownames(r$counts)]; nmat <- nmat[,rownames(r$counts)]; # restrict to cells that passed p2 filter
# take cluster labels
cluster.label <- r$clusters$PCA[[1]]
cell.colors <- pagoda2:::fac2col(cluster.label)
# take embedding
emb <- r$embeddings$PCA$tSNE

cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA)))

emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,
                                            cell.dist=cell.dist,
                                            fit.quantile=fit.quantile)

show.velocity.on.embedding.cor(emb,rvel.cd,n=300,scale='sqrt',
                               #cell.colors=ac(cell.colors,alpha=0.5),
                               cex=0.8,arrow.scale=5,show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,
                               do.par=F,cell.border.alpha = 0.1)

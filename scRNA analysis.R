library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggsci)
library(ClusterGVis)

logFCfilter=1               
adjPvalFilter=0.05          

#
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#
pbmc=CreateSeuratObject(counts = data,project = "seurat", min.cells=3, min.features=50, names.delim = "_")
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pbmc=subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 5)    

#
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)

#
pbmc=ScaleData(pbmc)         
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))    

#
pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:15)

#
pcSelect=15
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)       
pbmc <- FindClusters(object = pbmc, resolution = 0.5)        
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)           

#
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]

#
monocle.matrix=as.matrix(pbmc@assays$RNA@data)
monocle.sample=pbmc@meta.data
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.clusterAnn=clusterAnn
monocle.markers=sig.markers

#
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])

#
clusterAnn=as.character(monocle.clusterAnn[,2])
names(clusterAnn)=paste0("cluster",monocle.clusterAnn[,1])
pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),clusterAnn)

#
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, as.vector(sig.markers$gene))
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree')
cds <- orderCells(cds)

#
groups=subset(pData(cds),select='State')
pbmc=AddMetaData(object=pbmc, metadata=groups, col.name="group")
geneList=list()
for(i in levels(factor(groups$State))){
  pbmc.markers=FindMarkers(pbmc, ident.1 = i, group.by = 'group')
  sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
  sig.markers=cbind(Gene=row.names(sig.markers), sig.markers)
  geneList[[i]]=row.names(sig.markers)
}

#
BEAM_res <- BEAM(cds[unionGenes,],branch_point = 1, cores = 2)
BEAM_res <-BEAM_res[order(BEAM_res$qval),]
BEAM_res <-BEAM_res[,c("gene_short_name", "pval", "qval")]
BEAM_genes <- BEAM_res[order(BEAM_res$qval), "gene_short_name"][1:100]

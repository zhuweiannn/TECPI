library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#GSVA
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
gsvaResult=gsva(data, 
                geneSets, 
                min.sz=10, 
                max.sz=500, 
                verbose=TRUE,
                parallel.sz=1)
gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)

#
Subtype=read.table(SubtypeFile, header=T, sep="\t", check.names=F, row.names=1)
gsvaResult=t(gsvaResult)
sameSample=intersect(row.names(gsvaResult), row.names(Subtype))
gsvaResult=gsvaResult[sameSample,,drop=F]
Subtype=Subtype[sameSample,,drop=F]
gsvaSubtype=cbind(gsvaResult, Subtype)
Project=gsub("(.*?)\\_.*", "\\1", rownames(gsvaSubtype))
gsvaSubtype=cbind(gsvaSubtype, Project)

#
adj.P.Val.Filter=0.05
allType=as.vector(gsvaSubtype$Subtype)
comp=combn(levels(factor(allType)), 2)

#
for(i in 1:ncol(comp)){
  treat=gsvaSubtype[gsvaSubtype$Subtype==comp[2,i],]
  con=gsvaSubtype[gsvaSubtype$Subtype==comp[1,i],]
  data=rbind(con, treat)
  Type=as.vector(data$Subtype)
  ann=data[,c(ncol(data), (ncol(data)-1))]
  data=t(data[,-c((ncol(data)-1), ncol(data))])
  design=model.matrix(~0+factor(Type))
  colnames(design)=levels(factor(Type))
  fit=lmFit(data, design)
  contrast=paste0(comp[2,i], "-", comp[1,i])
  cont.matrix=makeContrasts(contrast, levels=design)
  fit2=contrasts.fit(fit, cont.matrix)
  fit2=eBayes(fit2)
  allDiff=topTable(fit2,adjust='fdr',number=200000)
  allDiffOut=rbind(id=colnames(allDiff),allDiff)
  write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
  diffSig=allDiff[with(allDiff, (abs(logFC)>0.1 & adj.P.Val < adj.P.Val.Filter )), ]
  diffSigOut=rbind(id=colnames(diffSig),diffSig)
  write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
  bioCol=c("#BC3C29FF", "#0072B5FF")
  ann_colors=list()
  m6aCluCol=bioCol[1:length(levels(factor(allType)))]
  names(m6aCluCol)=levels(factor(allType))
  ann_colors[["Subtype"]]=m6aCluCol[c(comp[1,i], comp[2,i])]
  termNum=50    
  diffTermName=as.vector(rownames(diffSig))
  diffLength=length(diffTermName)
  if(diffLength<termNum){termNum=diffLength}
  hmGene=diffTermName[1:termNum]
  hmExp=data[hmGene,]
  pdf(file=paste0(contrast,".heatmap.KEGG.pdf"), width=15, height=10)
  pheatmap(hmExp, 
           annotation=ann,
           annotation_colors = ann_colors,
           color = colorRampPalette(c("#0072B5FF","white","#BC3C29FF"))(100), #表示热图颜色,(100)表示100个等级
           cluster_cols =F,
           cluster_rows = T,
           show_colnames = F,
           gaps_col=as.vector(cumsum(table(Type))),
           scale="row",
           fontsize = 15,
           fontsize_col = 15,
           fontsize_row = 12)
  dev.off()
}


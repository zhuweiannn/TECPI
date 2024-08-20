library(GSVA)
library(limma)
library(GSEABase)
library(ggpubr)
library(reshape2)

rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]
geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())

ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}

ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
write.table(ssgseaOut,file="ssgseaOut.txt",sep="\t",quote=F,col.names=F)

#
scoreFile="ssgseaOut.txt"       
data=read.table(scoreFile,sep="\t",header=T,check.names=F,row.names=1)      

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(data))
data=avereps(t(data))

Subtype=read.table(SubtypeFile,header=T,sep="\t",row.names=1,check.names=F)

sameSample=intersect(row.names(data),row.names(Subtype))
data=data[sameSample,]
Subtype=Subtype[sameSample,]
rt=cbind(data,Subtype[,"Subtype"])
rt=rt[,-(ncol(rt)-1)]


immCell=c("aDCs","B_cells","CD8+_T_cells","DCs","iDCs","Macrophages",
          "Mast_cells","Neutrophils","NK_cells","pDCs","T_helper_cells",
          "Tfh","Th1_cells","Th2_cells","TIL","Treg")
rt1=rt[,c(immCell,"Subtype")]
data=melt(rt1,id.vars=c("Subtype"))
colnames(data)=c("Subtype","Type","Score")
data$Subtype=factor(data$Subtype, levels=c("C1","C2"))
boxplot=ggboxplot(data, x="Type", y="Score", fill = "Subtype",
                  orientation="horizontal",
                  xlab="",
                  ylab="Score",
                  legend.title="Subtype",
                  width=1,
                  font.x = c(18, "plain", "black"), 
                  font.y = c(18, "plain", "black"),  
                  palette=c("#BC3C29FF", "#0072B5FF"))+
  rotate_y_text(0)+
  stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif",size = 7, position = "identity", 
                     vjust = 0.7)+
  theme(legend.title=element_text(size=18),
        legend.text=element_text(size=15),
        axis.text = element_text(size=15),
        axis.title = element_text(size=15),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

#
immFunction=c("APC_co_inhibition","APC_co_stimulation","CCR",
              "Check-point","Cytolytic_activity","HLA","Inflammation-promoting",
              "MHC_class_I","Parainflammation","T_cell_co-inhibition",
              "T_cell_co-stimulation","Type_I_IFN_Reponse","Type_II_IFN_Reponse")
rt1=rt[,c(immFunction,"Subtype")]
data=melt(rt1,id.vars=c("Subtype"))
colnames(data)=c("Subtype","Type","Score")
data$Subtype=factor(data$Subtype, levels=c("C1","C2"))
boxplot=ggboxplot(data, x="Type", y="Score", fill = "Subtype",
                  orientation="horizontal",
                  xlab="",
                  ylab="Score",
                  legend.title="Subtype",
                  width=1,
                  font.x = c(18, "plain", "black"), 
                  font.y = c(18, "plain", "black"), 
                  palette=c("#BC3C29FF", "#0072B5FF"))+
  rotate_y_text(0)+
  stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif",size = 7, position = "identity", 
                     vjust = 0.7)+
  theme(legend.title=element_text(size=18),
        legend.text=element_text(size=15),
        axis.text = element_text(size=15),
        axis.title = element_text(size=15),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))



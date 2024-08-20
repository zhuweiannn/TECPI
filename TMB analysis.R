library(maftools)   
library(ggpubr)
library(reshape2)

#
TECPI=read.table("TECPI.txt", header=T, sep="\t", check.names=F)
outTab=TECPI[,c(1, ncol(TECPI))]
colnames(outTab)=c("Tumor_Sample_Barcode", "TECPI")

geneNum=20   
geneMut=read.table("geneMut.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneMut)[1:geneNum]

#
ann_colors=list()
col=c("#0072B5FF", "#BC3C29FF")
names(col)=c("Low", "High")
ann_colors[["TECPI"]]=col

# 
vc_cols = c("#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF", "#B09C85FF", "#3C5488FF")  
names(vc_cols) = c(
  'Missense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Del',
  'Nonsense_Mutation',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)


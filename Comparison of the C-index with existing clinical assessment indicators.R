library(dplyr)
library(survival)
library(rms)
library(pec)

TECPI=read.table(TECPIFile, header=T, sep="\t", check.names=F, row.names=1)
TECPI=TECPI[,c("BCR.time", "BCR", "TECPIScore")]
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

samSample=intersect(row.names(TECPI), row.names(cli))
TECPI1=TECPI[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(TECPI1, cli)
#
customCol <- c("#B24745FF", 
               "#374E55FF", 
               "#DF8F44FF", 
               "#00A1D5FF",  
               "#79AF97FF",  
               "#6A6599FF")  
bioCol <- customCol[1:(ncol(rt)-1)]

#
TECPIScore=cph(Surv(BCR.time,BCR)~TECPIScore, data=rt, surv=TRUE)
Age=cph(Surv(BCR.time,BCR)~Age, data=rt, surv=TRUE)
GleasonScore=cph(Surv(BCR.time,BCR)~GleasonScore, data=rt, surv=TRUE)
PSA=cph(Surv(BCR.time,BCR)~PSA, data=rt, surv=TRUE)
pTstage=cph(Surv(BCR.time,BCR)~pTstage, data=rt, surv=TRUE)
pNstage=cph(Surv(BCR.time,BCR)~pNstage, data=rt, surv=TRUE)
c_index  <- cindex(list("TECPI score"=TECPIScore, 
#                        "Age"=Age,
                        "Gleason Score"=GleasonScore,
                        "PSA"=PSA,
                        "pT stage"=pTstage),
                   formula=Surv(BCR.time,BCR)~ .,
                   data=rt,
                   eval.times=seq(0,15,1),
                   splitMethod="bootcv",
                   B=500
)

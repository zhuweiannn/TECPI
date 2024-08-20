library(limma)
library(ConsensusClusterPlus)

maxK=10
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=1000,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="pam",
                             distance="spearman")


#
clusterNum=2        #
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
cluster$cluster=paste0("C", cluster$cluster)
clusterOut=rbind(ID=colnames(cluster), cluster)



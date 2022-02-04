
#source("/mnt/ceph/jarredk/MRPC_UPDATE/pseudocode.R")
source("/mnt/ceph/jarredk/MRPC_UPDATE/MRPCgeneral.R")
library('MRPC', lib="/mnt/ceph/jarredk/Rpackages")
library('qvalue', lib="/mnt/ceph/jarredk/Rpackages")

XX=simu_data_layered


# PCs <- prcomp(XX,scale=TRUE)$x
# 
# All.Pvalues <- list()
# List.significant.asso1 <- list()
# 
# for (m in 1:(dim(PCs)[1]-1)) {
#   corr.PCs <- corr.test(PCs[,m],XX,use = 'pairwise.complete.obs', adjust = "none")
#   # The p values
#   Pvalues <- corr.PCs$p
#   Pvalues.nona <- Pvalues[!is.na(Pvalues)]
#   All.Pvalues [[m]] <- Pvalues.nona
#   qobj <- qvalue(Pvalues.nona, fdr.level=0.10) 
#   
#   # Significant associations
#   Significant.asso <- qobj$significant
#   List.significant.asso1[[m]] <- which(Significant.asso,useNames = TRUE)
# }
variant.type="eQTL"
V=XX[,1]
data=as.data.frame(XX[,-1])
U=NULL




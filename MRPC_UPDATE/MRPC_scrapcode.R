
#source("/mnt/ceph/jarredk/MRPC_UPDATE/pseudocode.R")
source("/mnt/ceph/jarredk/MRPC_UPDATE/MRPCgeneral.R")
library('MRPC', lib="/mnt/ceph/jarredk/Rpackages")
library('qvalue', lib="/mnt/ceph/jarredk/Rpackages")

XX=simu_data_layered

# top5=c(1,6, 33, 40, 48)
# path='/mnt/ceph/jarredk/Reg_Net/'
# tissues.vec=tissue.names[top5, 2:3]
# 
# t=1
# 
# #load additional files
# file1=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",
#             tissues.vec[t,1],"_AllPC/data.snp.cis.trans.final.",
#             tissues.vec[t,1],".V8.unique.snps.RData", sep = "")
# 
# file2=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",
#             tissues.vec[t,1],"_AllPC/PCs.matrix.",
#             tissues.vec[t,1], ".RData", sep="")
# 
# #load trio, pc, and confounders data
# trios=loadRData(fileName=file1)
# PCs=loadRData(fileName=file2)
# list.match.significant=loadRData("/mnt/ceph/jarredk/AddisReRunFiles/Updated_PC_Selection/List.Match.significant.trios.AdiposeSubcutaneous.RData")
# 
# #preallocate some test trios
# inputM0=cbind.data.frame(trios[,(2-1)*3+(1:3)], PCs[,list.match.significant[[2]]])
# inputM1.1=cbind.data.frame(trios[,(54-1)*3+(1:3)], PCs[,list.match.significant[[54]]])
# inputM1.2=cbind.data.frame(trios[,(510-1)*3+(1:3)], PCs[,list.match.significant[[510]]])
# inputM2=cbind.data.frame(trios[,(1214-1)*3+(1:3)], PCs[,list.match.significant[[1214]]])
# inputM3=cbind.data.frame(trios[,(1-1)*3+(1:3)], PCs[,list.match.significant[[1]]])
# inputM4=cbind.data.frame(trios[,(9-1)*3+(1:3)], PCs[,list.match.significant[[9]]])
# inputRARE=cbind.data.frame(trios[,(1071-1)*3+(1:3)], PCs[,list.match.significant[[1071]]])
# 
# #form a network with extra nodes
# data=cbind.data.frame(inputM1.1[,2:3], inputM2[,2:3])
# Z=cbind.data.frame(inputM1.1[,-c(1:3)],inputM2[,-c(1:3)])
# U=Z[match(unique(colnames(Z)), colnames(Z))]
# V=inputM1.1[,1]

PCs <- prcomp(XX,scale=TRUE)$x

All.Pvalues <- list()
List.significant.asso1 <- list()

for (m in 1:(dim(PCs)[1]-1)) {
  corr.PCs <- corr.test(PCs[,m],XX,use = 'pairwise.complete.obs', adjust = "none")
  # The p values
  Pvalues <- corr.PCs$p
  Pvalues.nona <- Pvalues[!is.na(Pvalues)]
  All.Pvalues [[m]] <- Pvalues.nona
  qobj <- qvalue(Pvalues.nona, fdr.level=0.10) 
  
  # Significant associations
  Significant.asso <- qobj$significant
  List.significant.asso1[[m]] <- which(Significant.asso,useNames = TRUE)
}
variant.type="eQTL"
V=XX[,1]
data=as.data.frame(XX[,-1])
U=NULL








sim=read.csv(file="/mnt/ceph/jarredk/GMACanalysis/master_tables/TRIOS_imbalanced_genotypes_WB_Updated_12_10_2021.csv", header = T)
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
addis.pcs=list()

pc.info=as.data.frame(matrix(0, nrow = length(sim$Trio.Num), ncol=3))
colnames(pc.info)=c("Trio.Num", "Num.pcs", "updated.per.var")

list.pcs=loadRData(fileName="/mnt/ceph/jarredk/AddisReRunFiles/Updated_PC_Selection/List.Match.significant.trios.WholeBlood.RData")

num.pcs=lapply(list.pcs, length)



tissue="WholeBlood"
#read in badsha peerdata.v8 file
Peerdata.V8=loadRData(paste0("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/PEER_Files_V8/Peerdata.",tissue,".V8.RData"))

#preform PCA in same method as badsha and retain eigen values
Peerdata2 <- t(Peerdata.V8)
Peerdata3 <- Peerdata2[,apply(Peerdata2, 2, var, na.rm=TRUE) != 0]
# pca matrix
PCs<- prcomp(Peerdata3,scale=TRUE)

evs=PCs$sdev^2
names(evs)=colnames(PCs$x)

for(i in 1:length(sim$Trio.Num)){
  
  pc.info$Trio.Num[i]=sim$Trio.Num[i]
  pc.info$Num.pcs[i]=num.pcs[[sim$Trio.Num[i]]]
  pc.info$updated.per.var[i]=sum(evs[match(colnames(PCs$x)[list.pcs[[sim$Trio.Num[i]]]], names(evs))])/sum(evs)
  
}


save(pc.info, "/mnt/ceph/jarredk/AddisReRunFiles/Updated_PC_Selection/Updated.PC.info.WholeBlood.RData")
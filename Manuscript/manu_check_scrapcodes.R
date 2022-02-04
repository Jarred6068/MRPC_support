
#load in trio master summary tables:
load("/mnt/ceph/jarredk/Manuscript/Tables_extra/All_Trios_Master_table_LOND.RData")
load("/mnt/ceph/jarredk/Manuscript/Tables_extra/All_Trios_Master_table_ADDIS.RData")
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")


#check the description of ICOSLG as mediator and summary statistics:

gen1=subset(ms2, cis.gene.name=="ICOSLG" & trans.gene=="ENSG00000277117")
indices=gen1[which(gen1$model=="M1"),]$Trio.Index
tissues=gen1[which(gen1$model=="M1"),]$Tissue

cors=matrix(0, nrow = length(indices), ncol = 3)
colnames(cors)=c("SNP:cis","SNP:trans","cis:trans")
for(i in 1:length(indices)){
  
  G1=Lond2Addis.lookup(trio.index=indices[i], tissue.name=tissues[i], with.pc=TRUE)
  ct=cor(G1$trio.with.pc[,1:3], use = "pairwise.complete.obs")
  cors[i,]=c(ct[1,2:3], ct[2,3])
}

print(median(abs(cors[,1])))
print(median(abs(cors[,2])))
print(median(abs(cors[,3])))


#check the description of MAP3K2-DT as mediator and summary statistics:

gen2=subset(ms2, cis.gene.name=="MAP3K2-DT" & model=="M1")
indices2=gen2$Trio.Index
tissues2=gen2$Tissue

cors2=matrix(0, nrow = length(indices2), ncol = 3)
colnames(cors2)=c("SNP:cis","SNP:trans","cis:trans")
for(i in 1:length(indices2)){
  
  G1=Lond2Addis.lookup(trio.index=indices2[i], tissue.name=tissues2[i], with.pc=TRUE)
  ct=cor(G1$trio.with.pc[,1:3], use = "pairwise.complete.obs")
  cors2[i,]=c(ct[1,2:3], ct[2,3])
}

print(median(abs(cors2[,1])))
print(median(abs(cors2[,2])))
print(median(abs(cors2[,3])))


#checking S11,S12,S4,S5
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")

M1.Data=loadRData(fileName = paste("/mnt/ceph/jarredk/Reg_Net/AL_M1Data/", "WholeBlood",".Rdata", sep = ""))
addisM1 = M1.Data$Catalog[[1]][[1]][[2]]
londM1 = M1.Data$Catalog[[1]][[1]][[1]]
Lond2Addis.lookup(trio.index=londM1$type2[3], tissue.name="WholeBlood", with.pc=TRUE)


load(file="/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/AdiposeSubcutaneous_AllPC/data.snp.cis.trans.final.AdiposeSubcutaneous.V8.unique.snps.RData")
#updated pc selection list
load("/mnt/ceph/jarredk/AddisReRunFiles/Updated_PC_Selection/List.significant.asso1.AdiposeSubcutaneous.RData")
load(file="/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/AdiposeSubcutaneous_AllPC/PCs.matrix.AdiposeSubcutaneous.RData")


# columns number for trios
start.col <- seq(1,dim(data.snp.cis.trans.final.AdiposeSubcutaneous.V8.unique.snps)[2],3)
end.col <- seq(3,dim(data.snp.cis.trans.final.AdiposeSubcutaneous.V8.unique.snps)[2],3)

trios <- dim(data.snp.cis.trans.final.AdiposeSubcutaneous.V8.unique.snps)[2]/3
print(trios)

Match.PCs.col <- list()
List.Match.significant.trios <- list()

for (ii in 1:trios) {
  # data
  data <- data.snp.cis.trans.final.AdiposeSubcutaneous.V8.unique.snps[,start.col[ii]:end.col[ii]]
  # search in each PCs
  for (m1 in 1:(dim(PCs.matrix.AdiposeSubcutaneous)[1]-1)) {
    Match.PCs.col[[m1]] <- intersect(colnames(data),colnames(data.snp.cis.trans.final.AdiposeSubcutaneous.V8.unique.snps)[List.significant.asso1.AdiposeSubcutaneous[[m1]]])
  }
  List.Match.significant <- which(sapply(Match.PCs.col, function(e) length(e)!=0))
  List.Match.significant.trios[[ii]] <- List.Match.significant
  if(length(List.Match.significant)!=0)
  {
    data.withPC <- cbind(data,PCs.matrix.AdiposeSubcutaneous[,List.Match.significant])
    colnames(data.withPC)[4:ncol(data.withPC)] <- paste("PC",List.Match.significant,sep="")
  }
  else
  {
    data.withPC <- data 
  }
}

save(List.Match.significant.trios, file="/mnt/ceph/jarredk/AddisReRunFiles/Updated_PC_Selection/List.Match.significant.trios.AdiposeSubcutaneous.RData")
  
  











#for WholeBlood


load(file="/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/data.snp.cis.trans.final.WholeBlood.V8.unique.snps.RData")
#updated pc selection list
List.significant.asso1.WholeBlood=loadRData("/mnt/ceph/jarredk/AddisReRunFiles/Updated_PC_Selection/List.significant.asso1.WholeBlood.RData")
load(file="/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/PCs.matrix.WholeBlood.RData")


# columns number for trios
start.col <- seq(1,dim(data.snp.cis.trans.final.WholeBlood.V8.unique.snps)[2],3)
end.col <- seq(3,dim(data.snp.cis.trans.final.WholeBlood.V8.unique.snps)[2],3)

trios <- dim(data.snp.cis.trans.final.WholeBlood.V8.unique.snps)[2]/3
print(trios)

Match.PCs.col <- list()
List.Match.significant.trios <- list()

for (ii in 1:trios) {
  # data
  data <- data.snp.cis.trans.final.WholeBlood.V8.unique.snps[,start.col[ii]:end.col[ii]]
  # search in each PCs
  for (m1 in 1:(dim(PCs.matrix.WholeBlood)[1]-1)) {
    Match.PCs.col[[m1]] <- intersect(colnames(data),colnames(data.snp.cis.trans.final.WholeBlood.V8.unique.snps)[List.significant.asso1.WholeBlood[[m1]]])
  }
  List.Match.significant <- which(sapply(Match.PCs.col, function(e) length(e)!=0))
  List.Match.significant.trios[[ii]] <- List.Match.significant
  if(length(List.Match.significant)!=0)
  {
    data.withPC <- cbind(data,PCs.matrix.WholeBlood[,List.Match.significant])
    colnames(data.withPC)[4:ncol(data.withPC)] <- paste("PC",List.Match.significant,sep="")
  }
  else
  {
    data.withPC <- data 
  }
}

save(List.Match.significant.trios, file="/mnt/ceph/jarredk/AddisReRunFiles/Updated_PC_Selection/List.Match.significant.trios.WholeBlood.RData")

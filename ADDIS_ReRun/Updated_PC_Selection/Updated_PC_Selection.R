

source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
#get the top five tissues by sample size
top5=c(1, 6, 33, 40, 48)
path='/mnt/ceph/jarredk/Reg_Net/'
tissues.vec=tissue.names[top5, 2]
library('psych', lib="/mnt/ceph/jarredk/Rpackages")

for(i in 1:length(tissues.vec)){

  print(paste0("Loading Data For...",tissues.vec[i]))
  PCs.matrix=loadRData(fileName=paste0("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissues.vec[i],"_AllPC/PCs.matrix.",tissues.vec[i],".RData"))
  data.all.unique.snps=loadRData(fileName=paste0("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissues.vec[i],
                                                 "_AllPC/data.snp.cis.trans.final.",tissues.vec[i],".V8.unique.snps.RData"))


  colnames(data.all.unique.snps) <- make.unique(colnames(data.all.unique.snps)) #duplicated colnames

  All.Pvalues <- list()
  List.significant.asso1 <- list()

  print(paste0("Getting PC List For...",tissues.vec[i]))
  # calculate associated PCs on all trios
  for (m in 1:(dim(PCs.matrix)[1]-1)) {
    corr.PCs <- corr.test(PCs.matrix[,m],data.all.unique.snps,use = 'pairwise.complete.obs', adjust="none")
    # The p values
    Pvalues <- corr.PCs$p
    Pvalues.nona <- Pvalues[!is.na(Pvalues)]
    All.Pvalues [[m]] <- Pvalues.nona
    qobj <- qvalue(Pvalues.nona, fdr.level=0.10)

    # Significant associations
    Significant.asso <- qobj$significant
    List.significant.asso1[[m]] <- which(Significant.asso,useNames = TRUE)
  }

  print(paste0("Saving List For...",tissues.vec[i]))
  save(List.significant.asso1, file=paste0("/mnt/ceph/jarredk/AddisReRunFiles/Updated_PC_Selection/List.significant.asso1.",tissues.vec[i],".RData"))
}




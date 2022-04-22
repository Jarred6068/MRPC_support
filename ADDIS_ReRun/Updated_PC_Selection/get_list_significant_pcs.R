

source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")

#get the top five tissues by sample size
top5=c(1, 6, 33, 40, 48)
path='/mnt/ceph/jarredk/Reg_Net/'
tissues.vec=tissue.names[, 2]

for(t in 1:length(tissues.vec)){

  PCs.matrix=loadRData(fileName=paste0("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissues.vec[t],"_AllPC/PCs.matrix.",tissues.vec[t],".RData"))
  data.all.unique.snps=loadRData(fileName=paste0("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissues.vec[t],
                                                 "_AllPC/data.snp.cis.trans.final.",tissues.vec[t],".V8.unique.snps.RData"))
  #updated pc selection list
  List.significant.asso1=loadRData(fileName=paste0("/mnt/ceph/jarredk/AddisReRunFiles/Updated_PC_Selection/List.significant.asso1.",tissues.vec[t],".RData"))


  # columns number for trios
  start.col <- seq(1,dim(data.all.unique.snps)[2],3)
  end.col <- seq(3,dim(data.all.unique.snps)[2],3)

  trios <- dim(data.all.unique.snps)[2]/3
  print(trios)

  Match.PCs.col <- list()
  List.Match.significant.trios <- list()
  data.sets=list()

  print(paste0("Matching For...",tissues.vec[t]))
  for (ii in 1:trios) {
    # data
    data <- data.all.unique.snps[,start.col[ii]:end.col[ii]]
    # search in each PCs
    for (m1 in 1:(dim(PCs.matrix)[1]-1)) {
      Match.PCs.col[[m1]] <- intersect(colnames(data),colnames(data.all.unique.snps)[List.significant.asso1[[m1]]])
    }
    List.Match.significant <- which(sapply(Match.PCs.col, function(e) length(e)!=0))
    List.Match.significant.trios[[ii]] <- List.Match.significant
    if(length(List.Match.significant)!=0)
    {
      data.withPC <- cbind(data,PCs.matrix[,List.Match.significant])
      colnames(data.withPC)[4:ncol(data.withPC)] <- paste("PC",List.Match.significant,sep="")
    }
    else
    {
      data.withPC <- data
    }
    data.sets[[ii]]=data.withPC
  }

  save(List.Match.significant.trios, file=paste0("/mnt/ceph/jarredk/AddisReRunFiles/Updated_PC_Selection/List.Match.significant.trios.",tissues.vec[t],".RData"))
  save(data.sets, file=paste0("/mnt/ceph/jarredk/AddisReRunFiles/Updated_PC_Selection/data.with.PCs.",tissues.vec[t],".RData"))

}











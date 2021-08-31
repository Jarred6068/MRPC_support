#GMAC analysis:

#import tables for top 5 tissues with largest sample sizes
#Adipose - Subcutaneous                     
#Artery - Tibial                            
#Skin - Sun Exposed (Lower leg)             
#Whole Blood                                
#Muscle - Skeletal             
source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")

top5=c(1,6, 40, 48, 33)
path='/mnt/ceph/jarredk/Reg_Net/'
tissues.vec=tissue.names[top5, 2:3]


run.GMAC=function(tissues.vec=tissue.vec, path.tables=path){
  
  types=c('LM1T1','LM1T2','AM1T1', 'LM1T2')

  for(t in 1:length(tissues.vec)){
  

    lm1t1=read.csv(paste0(path.tables, types[1], '/', tissues.vec[t,1], '.csv'), header = T)
    lm1t2=read.csv(paste0(path.tables, types[2], '/', tissues.vec[t,1], '.csv'), header = T)
    am1t1=read.csv(paste0(path.tables, types[3], '/', tissues.vec[t,1], '.csv'), header = T)
    am1t2=read.csv(paste0(path.tables, types[4], '/', tissues.vec[t,1], '.csv'), header = T)
  
  
    #load additional files
    file1=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",
                tissues.vec[t,1],"_AllPC/data.snp.cis.trans.final.",
                tissues.vec[t,1],".V8.unique.snps.RData", sep = "")
    
    file2=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",
                tissues.vec[t,1],"_AllPC/PCs.matrix.",
                tissues.vec[t,1], ".RData", sep="")
    
    # file3=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissues.vec[t],
    #             "_AllPC/List.significant.asso1.", tissues.vec[t,1], ".RData", sep="")
    # 
    # file4=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/PEER_Files_V8/Peerdata.", 
    #             tissues.vec[t,1], ".V8.RData", sep="")

  
    trios=loadRData(fileName=file1)
    PCs=loadRData(fileName=file2)
    #sig.asso.pcs=loadRData(fileName=file3)
    #edata=loadRData(fileName=file4)
    confounders=read.table(paste0("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/GTEx_Analysis_v8_eQTL_covariates/", 
                                 tissues.vec[t,2], '.v8.covariates.txt'), header = T, sep="\t", row.names = 1)[66:68,]

    print(paste('data loaded for', tissues.vec[t,2], sep = ' '))
    
    lm1t1.gmac.list=assemble.tables(trios, lm1t1, PCs, confounders)
    lm1t2.gmac.list=assemble.tables(trios, lm1t2, PCs, confounders)
    am1t1.gmac.list=assemble.tables(trios, am1t1, PCs, confounders)
    am1t2.gmac.list=assemble.tables(trios, am1t2, PCs, confounders)
    
    
    output.lm1t1 <- gmac(known.conf = lm1t1.gmac.list$known.conf, cov.pool = lm1t1.gmac.list$cov.pool, 
                         exp.dat = lm1t1.gmac.list$exp.dat, snp.dat.cis = lm1t1.gmac.list$snp.dat.cis, 
                         trios.idx = lm1t1.gmac.list$trios.idx, nperm = 10000, nominal.p = TRUE)
    print('lm1t1...done...saving')
    save(output.lm1t1, file = paste0('/mnt/ceph/jarredk/GMACanalysis/', tissues.vec[t,1], '/lm1t1.output.Rdata'))
    print('...done')
    
    output.lm1t2 <- gmac(known.conf = lm1t2.gmac.list$known.conf, cov.pool = lm1t2.gmac.list$cov.pool, 
                         exp.dat = lm1t2.gmac.list$exp.dat, snp.dat.cis = lm1t2.gmac.list$snp.dat.cis, 
                         trios.idx = lm1t2.gmac.list$trios.idx, nperm = 10000, nominal.p = TRUE)
    
    print('lm1t2...done...saving')
    save(output.lm1t2, file = paste0('/mnt/ceph/jarredk/GMACanalysis/', tissues.vec[t,1], '/lm1t2.output.Rdata'))
    print('...done')
    
    output.am1t1 <- gmac(known.conf = am1t1.gmac.list$known.conf, cov.pool = am1t1.gmac.list$cov.pool, 
                         exp.dat = am1t1.gmac.list$exp.dat, snp.dat.cis = am1t1.gmac.list$snp.dat.cis, 
                         trios.idx = am1t1.gmac.list$trios.idx, nperm = 10000, nominal.p = TRUE)
    
    print('am1t1...done...saving')
    save(output.am1t1, file = paste0('/mnt/ceph/jarredk/GMACanalysis/', tissues.vec[t,1], '/am1t1.output.Rdata'))
    print('...done')
    
    output.am1t2 <- gmac(known.conf = am1t2.gmac.list$known.conf, cov.pool = am1t2.gmac.list$cov.pool, 
                         exp.dat = am1t2.gmac.list$exp.dat, snp.dat.cis = am1t2.gmac.list$snp.dat.cis, 
                         trios.idx = am1t2.gmac.list$trios.idx, nperm = 10000, nominal.p = TRUE)
    
    print('am1t2...done...saving')
    save(output.am1t2, file = paste0('/mnt/ceph/jarredk/GMACanalysis/', tissues.vec[t,1], '/am1t2.output.Rdata'))
    print('...done')
    
    
    
    
  }
  
}



assemble.tables=function(trio.table, meta.table, PCmat, kc, seed=222){
  
  set.seed(seed)
  
  print("---SNPs---")
  SNPmat=as.data.frame(trio.table[,match(meta.table$SNP, colnames(trio.table))])
  print(dim(SNPmat))
  omit=attr(na.omit(SNPmat),"na.action")
  if (length(omit)>0){
    
    print(length(omit))
    SNPfactors=apply.cfactor(SNPmat)
    SNPmat=imputeMCA(SNPfactors, seed = seed)$completeObs
    SNPmat=as.data.frame(apply(SNPmat, 2, as.numeric))
    #print(SNPmat)
    
  }
  
  print("---Cis-Expression---")
  cis.exp=as.data.frame(trio.table[,match(meta.table$Cis.Gene.ID, colnames(trio.table))])
  print(dim(cis.exp))
  omit=attr(na.omit(cis.exp),"na.action")
  if (length(omit)>0){
    
    print(length(omit))
    imp=mice(cis.exp)
    cis.exp=complete(imp)
    print(cis.exp)
    
  }
  
  print("---Trans-Expression---")
  trans.exp=as.data.frame(trio.table[,match(meta.table$Trans.Gene.ID, colnames(trio.table))])
  print(dim(trans.exp))
  omit=attr(na.omit(trans.exp),"na.action")
  if (length(omit)>0){
    
    print(length(omit))
    imp=mice(trans.exp)
    trans.exp=complete(imp)
    print(trans.exp)
    
  }
  
  trio.idx=cbind.data.frame(c(1:dim(SNPmat)[2]), c(1:dim(cis.exp)[2]), c( (dim(cis.exp)[2]+1):(dim(cis.exp)[2]*2)) )
  colnames(trio.idx)=c('snp', 'cis', 'trans')
  
  
  gmac.list=list(kc, t(PCmat), t(cbind.data.frame(cis.exp, trans.exp)), t(SNPmat), trio.idx)
  names(gmac.list)=c("known.conf", "cov.pool", "exp.dat", "snp.dat.cis", "trios.idx")
  
  
  return(gmac.list)
  
}






apply.cfactor=function(dataframe){
  
  for(i in 1:dim(dataframe)[2]){
    
    dataframe[,i]=as.factor(dataframe[,i])
    
  }
  
  return(dataframe)
  
}




































































































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
  
  for(t in 3:length(tissues.vec[,1])){
    
    
    lm1t1=read.csv(paste0(path.tables, types[1], '/', tissues.vec[t,1], '.csv'), header = T)
    lm1t2=read.csv(paste0(path.tables, types[2], '/', tissues.vec[t,1], '.csv'), header = T)
    am1t1=read.csv(paste0(path.tables, types[3], '/', tissues.vec[t,1], '.csv'), header = T)
    am1t2=read.csv(paste0(path.tables, types[4], '/', tissues.vec[t,1], '.csv'), header = T)
    
    
    tables.c=rbind.data.frame(lm1t1, lm1t2, am1t1, am1t2)
    
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
    
    tables.gmac.list=assemble.tables(trios, tables.c, PCs, confounders)

    
    
    output <- gmac(known.conf = tables.gmac.list$known.conf, cov.pool = tables.gmac.list$cov.pool, 
                         exp.dat = tables.gmac.list$exp.dat, snp.dat.cis = tables.gmac.list$snp.dat.cis, 
                         trios.idx = tables.gmac.list$trios.idx, nperm = 10000, nominal.p = TRUE)
    
    out.table=cbind.data.frame(tables.c, output[[1]], output[[2]])
    colnames(out.table)=c(colnames(tables.c), 
                          paste0('pval_', colnames(output[[1]])), 
                          paste0('beta_change_', colnames(output[[2]])))
    
    out.list=list(out.table, tables.gmac.list, output[[3]])
    names(out.list)=c("output.table", "input.list", "cov.indicator.list")
    
    print('output...done...saving')
    save(out.list, file = paste0('/mnt/ceph/jarredk/GMACanalysis/', tissues.vec[t,1], '/output.Rdata'))
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
  
  
  gmac.list=list(kc, t(PCmat), t(cbind.data.frame(cis.exp, trans.exp)), t(SNPmat), trio.idx, trio.table)
  names(gmac.list)=c("known.conf", "cov.pool", "exp.dat", "snp.dat.cis", "trios.idx", 'trio.table')
  
  
  return(gmac.list)
  
}






apply.cfactor=function(dataframe){
  
  for(i in 1:dim(dataframe)[2]){
    
    dataframe[,i]=as.factor(dataframe[,i])
    
  }
  
  return(dataframe)
  
}




































































































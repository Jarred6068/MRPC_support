#GMAC analysis:

#import tables for top 5 tissues with largest sample sizes
#Adipose - Subcutaneous                     
#Artery - Tibial                            
#Skin - Sun Exposed (Lower leg)             
#Whole Blood                                
#Muscle - Skeletal             
source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")

top5=c(1,6, 33, 40, 48)
path='/mnt/ceph/jarredk/Reg_Net/'
tissues.vec=tissue.names[top5, 2:3]


run.GMAC=function(tissues.vec=tissue.vec, path.tables=path, mediation.type='cis', which.trios="all"){
  
  #types=c('LM1T1','LM1T2','AM1T1', 'LM1T2')
  
  for(t in 1:length(tissues.vec[,1])){
    
    #load additional files
    file1=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",
                tissues.vec[t,1],"_AllPC/data.snp.cis.trans.final.",
                tissues.vec[t,1],".V8.unique.snps.RData", sep = "")
    
    file2=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",
                tissues.vec[t,1],"_AllPC/PCs.matrix.",
                tissues.vec[t,1], ".RData", sep="")
    
    #load trio, pc, and confounders data
    trios=loadRData(fileName=file1)
    PCs=loadRData(fileName=file2)
    #sig.asso.pcs=loadRData(fileName=file3)
    #edata=loadRData(fileName=file4)
    confounders=read.table(paste0("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/GTEx_Analysis_v8_eQTL_covariates/", 
                                  tissues.vec[t,2], '.v8.covariates.txt'), header = T, sep="\t", row.names = 1)[66:68,]
    
    print(paste('data loaded for', tissues.vec[t,2], sep = ' '))
    
    #run assemble tables function
    tables.gmac.list=assemble.tables(trio.table=trios, PCmat = PCs, kc = confounders, which.imp = 'miss')

    print('data.assembled...running.GMAC')
    #run gmac analysis
    if(which.trios=="all"){
      
      if(mediation.type=='trans'){
        
        output <- gmac(known.conf = tables.gmac.list$known.conf, cov.pool = tables.gmac.list$cov.pool, 
                       exp.dat = tables.gmac.list$exp.dat, snp.dat.cis = tables.gmac.list$snp.dat.cis, 
                       trios.idx = tables.gmac.list$trios.idx[,c(1,3,2)], nperm = 10000, nominal.p = TRUE)
        
        
      }else{
        
        output <- gmac(known.conf = tables.gmac.list$known.conf, cov.pool = tables.gmac.list$cov.pool, 
                       exp.dat = tables.gmac.list$exp.dat, snp.dat.cis = tables.gmac.list$snp.dat.cis, 
                       trios.idx = tables.gmac.list$trios.idx, nperm = 10000, nominal.p = TRUE)
        
        
      }
      
    }else{
      
      
      if(mediation.type=='trans'){
        
        output <- gmac(known.conf = tables.gmac.list$known.conf, cov.pool = tables.gmac.list$cov.pool, 
                       exp.dat = tables.gmac.list$exp.dat, snp.dat.cis = tables.gmac.list$snp.dat.cis, 
                       trios.idx = tables.gmac.list$trios.idx[which.trios,c(1,3,2)], nperm = 10000, nominal.p = TRUE)
        
        
      }else{
        
        output <- gmac(known.conf = tables.gmac.list$known.conf, cov.pool = tables.gmac.list$cov.pool, 
                       exp.dat = tables.gmac.list$exp.dat, snp.dat.cis = tables.gmac.list$snp.dat.cis, 
                       trios.idx = tables.gmac.list$trios.idx[which.trios,], nperm = 10000, nominal.p = TRUE)
        
        
      }
      
    }
    
    
    
    #reorganize output for saving
    out.table=cbind.data.frame(tables.gmac.list$trio.ref[which.trios,], output[[1]], output[[2]])
    colnames(out.table)=c(colnames(tables.gmac.list$trio.ref), 
                          paste0('pval_', colnames(output[[1]])), 
                          paste0('effect_change_', colnames(output[[2]])))
    
    out.list=list(out.table, tables.gmac.list[1:5], output[[3]])
    names(out.list)=c("output.table", "input.list", "cov.indicator.list")
    
    

    
    #saving output
    print('output...done...saving')
    save(out.list, file = paste0('/mnt/ceph/jarredk/GMACanalysis/', tissues.vec[t,1], '/all_trios_output_', mediation.type, 'stability_check.Rdata'))
    print('...done')
    

  }
  
}



#--------function to assemble data for GMAC function---------

assemble.tables=function(trio.table, PCmat, kc, seed=222, which.imp='mice'){
  
  trio.ref=matrix(colnames(trio.table),nrow=length(colnames(trio.table))/3, ncol=3, byrow=T)
  print(trio.ref[1:10,])
  
  SNPs=trio.table[,c(1:(dim(trio.table)[2]/3))*3-2]
  
  SNP.unique=as.data.frame(SNPs[,match(unique(colnames(SNPs)), colnames(SNPs))])
  
  if(length(attr(na.omit(SNP.unique), 'na.action')) > 0 ){
    
    if(which.imp=='mice'){
      
      SNP.unique=complete(mice(SNP.unique, defaultMethod="polyreg"))
      print("...imputation..complete")
      
    }else{
      
      SNP.factors=apply.cfactor(SNP.unique)
      imp=imputeMCA(SNP.factors, seed=seed)
      SNP.unique=as.data.frame(apply(imp$completeObs,2,as.numeric))
      #colnames(SNP.unique)=colnames(SNP.factors)
      print("...imputation..complete")
      print(SNP.unique[1:5,1:5])
      
    }
    
    
  }
  
  
  expression.mat = trio.table[, -(c(1:(dim(trio.table)[2]/3))*3-2)]
  express.unique = expression.mat[, match(unique(colnames(expression.mat)), colnames(expression.mat))]
  
  # if(length(attr(na.omit(express.unique), 'na.action')) > 0 ){
  #   
  #   express.unique=mice(express.unique, defaultMethod='pmm')
  #   
  # }
  
  
  trio.idx=build.triomap(trio.name.list = trio.ref, 
                         express = express.unique, 
                         genotype = SNP.unique)
  
  colnames(trio.idx)=c('snp', 'cis', 'trans')
  colnames(trio.ref)=c('snp', 'cis', 'trans')
  
  #assemble into list and transpose necessary tables so samples are in columns 
  gmac.list=list(kc, t(PCmat), t(express.unique), t(SNP.unique), trio.idx, trio.ref)
  names(gmac.list)=c("known.conf", "cov.pool", "exp.dat", "snp.dat.cis", "trios.idx", 'trio.ref')
  
  
  return(gmac.list)
  
}





# converts all columns of a matrix to factors 
apply.cfactor=function(dataframe){
  
  cnames=colnames(dataframe)
  
  for(i in 1:dim(dataframe)[2]){
    
    dataframe[,i]=as.factor(dataframe[,i])
    
  }
  
  colnames(dataframe)=cnames
  
  return(dataframe)
  
}


#a function to construct the trio index map in the expression and genotype dataframes

build.triomap = function(trio.name.list=NA, express, genotype){
  
  index.mat=as.data.frame(matrix(0, nrow=dim(trio.name.list)[1], ncol = 3))
  
  for(i in 1:dim(trio.name.list)[1]){
    
    
    index.mat[i,1] = match(trio.name.list[i,1], colnames(genotype))
    index.mat[i,2] = match(trio.name.list[i,2], colnames(express))
    index.mat[i,3] = match(trio.name.list[i,3], colnames(express))
    
  }
  
  return(index.mat)
  
}




































































































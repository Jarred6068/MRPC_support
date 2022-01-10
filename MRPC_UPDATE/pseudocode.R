

source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")
library('ppcor', lib="/mnt/ceph/jarredk/Rpackages")

top5=c(1,6, 33, 40, 48)
path='/mnt/ceph/jarredk/Reg_Net/'
tissues.vec=tissue.names[top5, 2:3]

t=1

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
list.match.significant=loadRData("/mnt/ceph/jarredk/AddisReRunFiles/Updated_PC_Selection/List.Match.significant.trios.AdiposeSubcutaneous.RData")

#preallocate some test trios
inputM0=cbind.data.frame(trios[,(2-1)*3+(1:3)], PCs[,list.match.significant[[2]]])
inputM1.1=cbind.data.frame(trios[,(54-1)*3+(1:3)], PCs[,list.match.significant[[54]]])
inputM1.2=cbind.data.frame(trios[,(510-1)*3+(1:3)], PCs[,list.match.significant[[510]]])
inputM2=cbind.data.frame(trios[,(1214-1)*3+(1:3)], PCs[,list.match.significant[[1214]]])
inputM3=cbind.data.frame(trios[,(1-1)*3+(1:3)], PCs[,list.match.significant[[1]]])
inputM4=cbind.data.frame(trios[,(9-1)*3+(1:3)], PCs[,list.match.significant[[9]]])
inputRARE=cbind.data.frame(trios[,(1071-1)*3+(1:3)], PCs[,list.match.significant[[1071]]])


#codd

MRPC.updated=function(data = NULL, nperms=1000, graph=TRUE){
  
  library('psych', "/mnt/ceph/jarredk/Rpackages")
  
  #allocate space for output
  output=list(marginal.tests=NULL, models=list(cis.model=NULL, trans.model=NULL), adj.matrix=NULL, 
              permutation=list(cis=NULL, trans=NULL), graph=NULL)
  
  #allocate adjacency matrix
  adj.mat=diag(0, nrow = 3, ncol = 3)
  colnames(adj.mat)=colnames(data)[1:3]
  row.names(adj.mat)=colnames(data)[1:3]
  
  
  #preform marginal tests with the SNP
  cor.result = corr.test(data[,1:3], use = "pairwise.complete.obs")
  output$marginal.tests=cor.result
  
  adj.mat[1,2:3]=replace(adj.mat[1,2:3], which(cor.result$p[1,2:3]<0.05), 1)
  #update adjacency matrix in output list
  output$adj.matrix=adj.mat
  
  #get the minor allele frequency
  minor.allele.freq=get.allele.freq(SNP=data[,1])
  
  if(minor.allele.freq>0.05){
    
    #preform standard regressions with minor allele freq. > 5%
    output=check.mediation(data=data, adj.mat = adj.mat, output = output )
    
  }else{
    
    print(paste0("Minor Allele Frequency < 5%...Preforming Permutation"))
    print(paste0("MAF:",minor.allele.freq))
    
    if(sum(adj.mat[1,2:3])==2){
      
      #preform permutation test with cis gene as mediator
      perm.result.cis=permuted.reg(data = data, nperms = nperms, med.type = "cis")
      #store result to adjacency
      adj.mat[2, 3]=replace(adj.mat[2,3], perm.result.cis$nominal.p<0.05, 1)
      #preform permutation test with trans gene as mediator
      perm.result.trans=permuted.reg(data = data, nperms = nperms, med.type = "trans")
      #store result to adjacency
      adj.mat[3, 2]=replace(adj.mat[3,2], perm.result.trans$nominal.p<0.05, 1)
      
      #store all results
      output$permutation$cis=perm.result.cis
      output$permutation$trans=perm.result.trans
      output$adj.matrix=adj.mat
      
    }else if (sum(adj.mat[1,2:3])==1){
      
      if(adj.mat[1,3]==1){
        
        #preform permutation test with trans gene as mediator
        perm.result.trans=permuted.reg(data = data, nperms = nperms, med.type = "trans")
        #store result to adjacency
        adj.mat[3, 2]=replace(adj.mat[3,2], perm.result.trans$nominal.p<0.05, 1)
        
        #store all results
        output$permutation$trans=perm.result.trans
        output$adj.matrix=adj.mat
        
      }else{
        
        #preform permutation test with cis gene as mediator
        perm.result.cis=permuted.reg(data = data, nperms = nperms, med.type = "cis")
        #store result to adjacency
        adj.mat[2, 3]=replace(adj.mat[2,3], perm.result.cis$nominal.p<0.05, 1)
        
        #store all results
        output$permutation$cis=perm.result.cis
        output$adj.matrix=adj.mat
      }
    }
    
    
  }
  
  
  print(output$adj.matrix)
  
  return(output)
  
}




#-----------------------Helper-function----------------------

get.allele.freq=function(SNP=NULL){
  
  #remove missing values
  SNP=na.omit(SNP)

  
  if(length(unique(SNP))>2){
    
    #get the genotype counts
    totals=summary(as.factor(SNP))
    
    #calculate allele frequencies
    reference=(2*totals[1]+totals[2])/(2*sum(totals))
    alternative=1-reference
    
  }else{
    #for SNPs with <3 genotypes
    
    if( length(which(names(summary(factor(SNP)))=="0"))==0 ){
      
      #calculate allele frequencies
      totals=c(NA ,summary(factor(SNP)))
      reference=totals[2]/(2*sum(totals))
      alternative=1-reference
      
    }else if ( length(which(names(summary(factor(SNP)))=="1"))==0 ){
      
      #calculate allele frequencies
      totals=c(summary(factor(SNP))[1], NA , summary(factor(SNP))[2])
      reference=(2*totals[1])/(2*sum(totals))
      alternative=1-reference
      
    }else{
      
      #calculate allele frequencies
      totals=c(summary(factor(SNP)), NA)
      reference=(2*totals[1]+totals[2])/(2*sum(totals))
      alternative=1-reference
    }
    
    
  }
  
  #return minor allele frequency
  return(minor.allele.freq=min(reference,alternative))
  
  
  
}




#-----------------------Helper-function-MRPC.updated----------------------

check.mediation=function(data=NULL, adj.mat=NULL, output=NULL){
  
  if(sum(adj.mat[1,2:3])==2){
    
    #preform regression with cis.gene as mediator
    cis.model=lm(data[,3]~., data=data[,-3])
    #store regression result
    cis.med.result=as.data.frame(summary(cis.model)$coefficients)
    #change adj. matrix based on test result
    adj.mat[2, 3]=replace(adj.mat[2,3], cis.med.result$`Pr(>|t|)`[3]<0.05, 1)
    
    #preform regression with trans gene as mediator
    trans.model=lm(data[,2]~., data=data[,-2])
    #store regression result
    trans.med.result=as.data.frame(summary(trans.model)$coefficients)
    #change adj. matrix based on test result
    adj.mat[3, 2]=replace(adj.mat[3,2], trans.med.result$`Pr(>|t|)`[3]<0.05, 1)
    
    output$models$cis.model=cis.model
    output$models$trans.model=trans.model
    output$adj.matrix=adj.mat
    
    
    
  }else if (sum(adj.mat[1,2:3])==1){
    
    if(adj.mat[1,3]==1){
      
      #preform regression with trans gene as mediator
      trans.model=lm(data[,2]~., data=data[,-2])
      #store regression result
      trans.med.result=as.data.frame(summary(trans.model)$coefficients)
      #change adj. matrix based on test result
      adj.mat[3, 2]=replace(adj.mat[3,2], trans.med.result$`Pr(>|t|)`[3]<0.05, 1)
      
      output$models$trans.model=trans.model
      output$adj.matrix=adj.mat
      
    }else{
      
      #preform regression with cis.gene as mediator
      cis.model=lm(data[,3]~., data=data[,-3])
      #store regression result
      cis.med.result=as.data.frame(summary(cis.model)$coefficients)
      #change adj. matrix based on test result
      adj.mat[2, 3]=replace(adj.mat[2,3], cis.med.result$`Pr(>|t|)`[3]<0.05, 1)
      
      output$models$cis.model=cis.model
      output$adj.matrix=adj.mat
      
    }
  }
  
  return(output)
  
}



#-----------------------Helper-function-for-permuted.reg----------------------


help.fn1=function(perm.map=NULL, data=NULL, med.type=NULL){
  
  #permute
  if(med.type=="cis"){
    
    data[,2]=data[,2][perm.map]
    
  }else{
    
    data[,3]=data[,3][perm.map]
    
  }
  
  
  #run regression
  if(med.type=="cis"){
    
    coef.mat=as.data.frame(summary(lm(data[,3]~., data=data[,-3]))$coefficients)
    wald.stat=coef.mat$`t value`[which(row.names(coef.mat)==colnames(data)[2])]
    #print(coef.mat)
    
  }else{
    
    coef.mat=as.data.frame(summary(lm(data[,2]~., data=data[,-2]))$coefficients)
    wald.stat=coef.mat$`t value`[which(row.names(coef.mat)==colnames(data)[3])]
    #print(coef.mat)
    
  }
  
  #return the wald stat
  return(wald.stat)
  
}





#-----------------------Permuted-Regression-function----------------------

permuted.reg=function(data=NULL, nperms=1000, med.type=NULL){
  
  #allocation of space
  wald.stat=NULL
  mediator_perm=matrix(c(1:dim(data)[1]), nrow=dim(data)[1], ncol = nperms)
  
  #get obs. wald statistic
  if(med.type=="cis"){
    
    coef.mat=as.data.frame(summary(lm(data[,3]~., data=data[,-3]))$coefficients)
    test.wald=coef.mat$`t value`[which(row.names(coef.mat)==colnames(data)[2])]
    
  }else{
    
    coef.mat=as.data.frame(summary(lm(data[,2]~., data=data[,-2]))$coefficients)
    test.wald=coef.mat$`t value`[which(row.names(coef.mat)==colnames(data)[3])]
    
  }
  
  
  #preallocate all permutations
  for (j in 0:2) {
    ind <- which(data[,1] == j)
    if (length(ind) > 1) {
      mediator_perm[ind, ] <- apply(mediator_perm[ind, ], 2, sample)
    }
  }

  #apply permutation regression
  wald.stat=apply(mediator_perm, 2, help.fn1, data=data, med.type=med.type)
  

  p.value=2*min(sum(wald.stat<test.wald)/length(wald.stat),sum(wald.stat>test.wald)/length(wald.stat) )
  
  nominal.p.value=2 * (1 - pnorm(abs((test.wald - mean(wald.stat))/sd(wald.stat))))
  
  
  return(list(p.value=p.value, nominal.p=nominal.p.value, null.wald.stats=wald.stat, obs.wald.stat=test.wald))
  
}



















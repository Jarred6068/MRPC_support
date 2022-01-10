
source("/mnt/ceph/jarredk/MRPC_UPDATE/pseudocode.R")


MRPC.updated2=function(data = NULL, inst.var=1, nperms=1000, graph=TRUE){
  
  library('psych', "/mnt/ceph/jarredk/Rpackages")
  
  #allocate space for output
  output=list(marginal.tests=NULL, models=list(cis.model=NULL, trans.model=NULL), adj.matrix=NULL, 
              permutation=list(cis=NULL, trans=NULL), graph=NULL)
  
  test.cases=list(matrix(0,1,0,0, nrow=2, ncol=2, byrow=T),matrix(0,0,0,1, nrow=2, ncol=2, byrow=T),
                  matrix(1,1,1,0, nrow=2, ncol=2, byrow=T),matrix(1,0,1,1, nrow=2, ncol=2, byrow=T),
                  matrix(0,1,0,1, nrow=2, ncol=2, byrow=T),matrix(1,1,1,1, nrow=2, ncol=2, byrow=T))
  
  #allocate adjacency matrix
  adj.mat=diag(0, nrow = 3, ncol = 3)
  colnames(adj.mat)=colnames(data)[1:3]
  row.names(adj.mat)=colnames(data)[1:3]
  
  #get the minor allele frequency
  minor.allele.freq=get.allele.freq(SNP=data[,1])
  
  if(minor.allele.freq>0.05){
    
    #preform standard regressions with minor allele freq. > 5%
    output=confirm.model(data=data, cases=test.cases, output = output )
    
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

confirm.model=function(data=NULL, adj.mat=NULL, output=NULL){
  
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



















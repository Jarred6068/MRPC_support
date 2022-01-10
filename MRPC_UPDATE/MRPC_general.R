
source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")
library('ppcor', lib="/mnt/ceph/jarredk/Rpackages")
library('psych', lib="/mnt/ceph/jarredk/Rpackages")

#-----------------------Helper-function-To-Get-Variant-Frequency---------------------

get.freq=function(V=NULL, type=c("eQTL", "CNV", "Other")){
  
  #remove missing values
  V=na.omit(V)
  
  if(type="eQTL"){
    if(length(unique(V))>2){
      
      #get the genotype counts
      totals=summary(as.factor(V))
      #calculate allele frequencies
      reference=(2*totals[1]+totals[2])/(2*sum(totals))
      alternative=1-reference
      
    }else{
      #for Vs with <3 genotypes
      
      if( length(which(names(summary(factor(V)))=="0"))==0 ){
        
        #calculate allele frequencies
        totals=c(NA ,summary(factor(V)))
        reference=totals[2]/(2*sum(totals))
        alternative=1-reference
        
      }else if ( length(which(names(summary(factor(V)))=="1"))==0 ){
        
        #calculate allele frequencies
        totals=c(summary(factor(V))[1], NA , summary(factor(V))[2])
        reference=(2*totals[1])/(2*sum(totals))
        alternative=1-reference
        
      }else{
        
        #calculate allele frequencies
        totals=c(summary(factor(V)), NA)
        reference=(2*totals[1]+totals[2])/(2*sum(totals))
        alternative=1-reference
      }
    }
    
    #return minor allele frequency
    return(minor.freq=min(reference,alternative))
      
  }else{
      
    freq=summary(as.factor(V))/sum(summary(as.factor(V)))
    #return minor allele frequency
    return(minor.freq=min(freq))
    
  }

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

permuted.reg=function(data=NULL, nperms=1000, med.type=c("T1","T2")){
  
  #allocation of space
  wald.stat=NULL
  mediator_perm=matrix(c(1:dim(data)[1]), nrow=dim(data)[1], ncol = nperms)
  
  #get obs. wald statistic
  if(med.type=="T1"){
    
    coef.mat=as.data.frame(summary(lm(data[,3]~., data=data[,-3]))$coefficients)
    test.wald=coef.mat$`t value`[which(row.names(coef.mat)==colnames(data)[2])]
    var.p=coef.mat$`Pr(>|t|)`[which(row.names(coef.mat)==colnames(data)[1])]
    
  }else{
    
    coef.mat=as.data.frame(summary(lm(data[,2]~., data=data[,-2]))$coefficients)
    test.wald=coef.mat$`t value`[which(row.names(coef.mat)==colnames(data)[3])]
    var.p=coef.mat$`Pr(>|t|)`[which(row.names(coef.mat)==colnames(data)[1])]
    
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
  
  
  #p.value=2*min(sum(wald.stat<test.wald)/length(wald.stat),sum(wald.stat>test.wald)/length(wald.stat) )
  
  nominal.p.value=2 * (1 - pnorm(abs((test.wald - mean(wald.stat))/sd(wald.stat))))
  
  
  return(list(nominal.p=nominal.p.value, variant.p=var.p))
  
}


#-------------------------Helper-fn-to-get-all-combinations-of-trios---------------------

allocate.all.trios=function(data=NULL, combs=NULL, with.variant=TRUE, V=NULL){
  
  #preallocate list
  trio.list=vector('list', length = length(combs))
  
  if(with.variant==TRUE){
    
    for(i in 1:length(trio.list)){
      
      #create n x 3 dataframes for each trio formed with the variant
      trio.list[[i]] = cbind.data.frame(V, data[,combs[1,i]], data[,combs[2,i]])
      
    }
    
  }else{
    
    for(i in 1:length(trio.list)){
      
      #create n x 3 dfs for each trio of nodes only
      trio.list[[i]] = cbind.data.frame(data[,combs[1,i]], data[,combs[2,i]], data[,combs[3,i]])
      
    }
    
  }

  
  return(trio.list)
  
}




#--------------------------Function-to-get-graph-skeleton--------------------------------

get.skel=function(V=NULL, data=NULL, U=NULL, p=NULL){
  
  #merge all variables
  X=cbind.data.frame(V,data,U)
  #get partial correlations obj
  H=pcor(X)$p.value
  #acquire adjacency based on pcor.test pvalues
  A=replace(H[,1:p], H[,1:p]<0.05, 1)
  
  return(A)
  
}


#------------------------wrapper-function---------------------------
MRPC.update=function(V=NULL, data=NULL, U=NULL, gamma=0.05, m=1000, variant.type=c("eQTL", "CNV", "Other")){
  
  #initialization
  p=dim(data)[2]
  n=dim(data)[1]
  q=dim(U)[2]
  
  combs.with.variant=combn(p, 2)
  combs.Tnodes=combn(p,3)
  
  # A=matrix(1, nrow = p+1, ncol = p+1)
  # diag(A)=rep(0, p)
  
  
  #get all trios
  trio.list.with.variant=allocate.all.trios(data=data, 
                                            combs=combs.with.variant, 
                                            with.variant=TRUE, 
                                            V=V)
  
  trio.list.with.variant=allocate.all.trios(data=data, 
                                            combs=combs.Tnodes, 
                                            with.variant=TRUE, 
                                            V=V)
  
  #get the graph skeleton
  A=get.skel(V=V, data=data, U=U)
  #calculate minor variant frequency
  variant.freq=get.freq(V=V, type=variant.type)
  
  if(variant.freq<gamma){
    
    out.coeffs.T1=lapply(trio.list.with.variant, permuted.reg, nperms=m, med.type="T1")
    out.coeffs.T2=lapply(trio.list.with.variant, permuted.reg, nperms=m, med.type="T2")
    
  }else{
    
    
  }
  
  
}



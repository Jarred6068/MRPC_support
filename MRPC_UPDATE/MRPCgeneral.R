
source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")
library('ppcor', lib="/mnt/ceph/jarredk/Rpackages")
library('psych', lib="/mnt/ceph/jarredk/Rpackages")




#-----------------------Helper-function-To-Get-Variant-Frequency---------------------

get.freq=function(V=NULL, type=c("eQTL", "CNV", "Other")){
  
  #remove missing values
  V=na.omit(V)
  
  if(type=="eQTL"){
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
  if(med.type=="T1"){
    
    data[,2]=data[,2][perm.map]
    
  }else{
    
    data[,3]=data[,3][perm.map]
    
  }
  
  
  #run regression
  if(med.type=="T1"){
    
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


#------------------------Helper-fn2-Regression----------------------------

do.reg=function(data=NULL, med.type=c("T1","T2")){
  
    #get obs. wald statistic
  if(med.type=="T1"){
      
    coef.mat=as.data.frame(summary(lm(data[,3]~., data=data[,-3]))$coefficients)
    test.wald=coef.mat$`t value`[which(row.names(coef.mat)==colnames(data)[2])]
    wald.p=coef.mat$`Pr(>|t|)`[which(row.names(coef.mat)==colnames(data)[2])]
    var.p=coef.mat$`Pr(>|t|)`[which(row.names(coef.mat)==colnames(data)[1])]
      
  }else{
      
    coef.mat=as.data.frame(summary(lm(data[,2]~., data=data[,-2]))$coefficients)
    test.wald=coef.mat$`t value`[which(row.names(coef.mat)==colnames(data)[3])]
    wald.p=coef.mat$`Pr(>|t|)`[which(row.names(coef.mat)==colnames(data)[3])]
    var.p=coef.mat$`Pr(>|t|)`[which(row.names(coef.mat)==colnames(data)[1])]
      
  }
  
  
  return(list(coef.mat=coef.mat, test.wald=test.wald, wald.p=wald.p, var.p=var.p))
  
}


#-----------------------Permuted-Regression-function----------------------

#runs both T1 as mediator and T2 as mediator and reports cor(V,T2)
reg.with.variant=function(data=NULL, permuted=TRUE, nperms=1000, return.indicator=FALSE, alpha=0.05){
  
  #allocation of space
  result=NULL
  wald.stat=NULL
  mediator_perm=matrix(c(1:dim(data)[1]), nrow=dim(data)[1], ncol = nperms)

  
  if(permuted==FALSE){
    
    #preform standard regressions
    T1.med=do.reg(data=data, med.type = "T1")
    T2.med=do.reg(data=data, med.type = "T2")
    
    ct=corr.test(data[,1:3], use="pairwise.complete.obs")  #set adjust = "none" to remove bonf.adj.
    result=c(T1.med$wald.p, T1.med$var.p, T2.med$wald.p, T2.med$var.p, ct$p[1,3])
    names(result)=c("B22", "B12", "B21", "B11", "V indep T2")
    
    #returns indicator values for significant tests
    if(return.indicator==TRUE){
      result=replace(result, result<alpha, 1)
      result=replace(result, result!=1, 0)
    }
    
    return(list(X=result, trio.name=colnames(data[,1:3])))
    
  }else{
      
    T1.med=do.reg(data=data, med.type = "T1")
    T2.med=do.reg(data=data, med.type = "T2")
      
    #preallocate all permutations
    for (j in 0:2) {
      ind <- which(data[,1] == j)
      if (length(ind) > 1) {
        mediator_perm[ind, ] <- apply(mediator_perm[ind, ], 2, sample)
      }
    }
      
    #apply permutation regression
    wald.stat.T1=apply(mediator_perm, 2, help.fn1, data=data, med.type="T1")
    wald.stat.T2=apply(mediator_perm, 2, help.fn1, data=data, med.type="T2")
      
      
    #p.value=2*min(sum(wald.stat<test.wald)/length(wald.stat),sum(wald.stat>test.wald)/length(wald.stat) )
    nominal.p.value.T1=2 * (1 - pnorm(abs((T1.med$test.wald - mean(wald.stat.T1))/sd(wald.stat.T2))))
    nominal.p.value.T2=2 * (1 - pnorm(abs((T2.med$test.wald - mean(wald.stat.T2))/sd(wald.stat.T2))))
    
    ct=corr.test(data[,1:3], use="pairwise.complete.obs")  #set adjust = "none" to remove bonf.adj.
    result=c(nominal.p.value.T1, T1.med$var.p, nominal.p.value.T2, T2.med$var.p, ct$p[1,3])
    names(result)=c("B22", "B12", "B21", "B11", "V indep T2")
    
    #returns indicator values for significant tests
    if(return.indicator==TRUE){
      result=replace(result, result<alpha, 1)
      result=replace(result, result!=1, 0)
    }
    
    return(list(X=result, trio.name=colnames(data[,1:3])))
    
  }
  
}


#-------------------------Helper-fn-to-get-all-combinations-of-trios---------------------

allocate.all.trios=function(data=NULL, combs=NULL, with.variant=TRUE, V=NULL, U=NULL){
  
  #preallocate list
  trio.list=vector('list', length = dim(combs)[2])
  
  if(with.variant==TRUE){
    
    for(i in 1:length(trio.list)){
      
      #create n x 3 dataframes for each trio formed with the variant
      if(is.null(U)){
        trio.list[[i]] = cbind.data.frame(V, data[,combs[1,i]], data[,combs[2,i]])
        colnames(trio.list[[i]])=c("V", names(data[,combs[1:2,i]]))
      }else{
        trio.list[[i]] = cbind.data.frame(V, data[,combs[1,i]], data[,combs[2,i]], U)
        colnames(trio.list[[i]])=c("V", names(data[,combs[1:2,i]]), colnames(U))
      }
      
    }
    
  }else{
    
    for(i in 1:length(trio.list)){
      
      #create n x 3 dfs for each trio of nodes only
      if(is.null(U)){
        trio.list[[i]] = cbind.data.frame(data[,combs[1,i]], data[,combs[2,i]], data[,combs[3,i]])
        colnames(trio.list[[i]])=c(colnames(data[,combs[1:3,i]]))
      }else{
        trio.list[[i]] = cbind.data.frame(data[,combs[1,i]], data[,combs[2,i]], data[,combs[3,i]], U)
        colnames(trio.list[[i]])=c(colnames(data[,combs[1:3,i]]), colnames(U))
      }
      
    }
    
  }

  
  return(trio.list)
  
}




#--------------------------Function-to-get-graph-skeleton--------------------------------

get.skel=function(data=NULL, p=NULL){
  
  if(sum(is.na(data))>0){
   X=na.omit(data)
   message("NA's detected in get.skel()...removing")
  }
  #get partial correlations obj
  H=pcor(X)$p.value
  #acquire adjacency based on pcor.test pvalues
  A=replace(H, H<0.05, 1)[1:p,1:p]
  B=replace(A, A!=1, 0)
  diag(B)=0
  
  return(B)
  
}



#--------------------------Helper-fn-for-reg.without.variant()---------------------
do.reg.wo.var=function(data=NULL, response=c("T1","T2","T3"), return.indicator=FALSE, alpha=0.05){

  
  #get obs. wald statistic
  if(response=="T1"){
    
    coef.mat=as.data.frame(summary(lm(data[,1]~., data=data[,-1]))$coefficients)
    #test.wald=coef.mat$`t value`[which(row.names(coef.mat)==colnames(data)[2])]
    wald.p.T2=coef.mat$`Pr(>|t|)`[which(row.names(coef.mat)==colnames(data)[2])]
    wald.p.T3=coef.mat$`Pr(>|t|)`[which(row.names(coef.mat)==colnames(data)[3])]
    
    wald.stats=c(wald.p.T2, wald.p.T3)
    names(wald.stats)=c("wald.p.T2","wald.p.T3")
    #returns indicator values for significant tests
    if(return.indicator==TRUE){
      wald.stats=replace(wald.stats, wald.stats<alpha, 1)
      wald.stats=replace(wald.stats, wald.stats!=1, 0)
    }
    
    return(list(coef.mat=coef.mat, pvalues=wald.stats))
    
  }else if(response=="T2"){
    
    coef.mat=as.data.frame(summary(lm(data[,2]~., data=data[,-2]))$coefficients)
    #test.wald=coef.mat$`t value`[which(row.names(coef.mat)==colnames(data)[3])]
    wald.p.T1=coef.mat$`Pr(>|t|)`[which(row.names(coef.mat)==colnames(data)[1])]
    wald.p.T3=coef.mat$`Pr(>|t|)`[which(row.names(coef.mat)==colnames(data)[3])]
    
    wald.stats=c(wald.p.T1, wald.p.T3)
    names(wald.stats)=c("wald.p.T1","wald.p.T3")
    #returns indicator values for significant tests
    if(return.indicator==TRUE){
      wald.stats=replace(wald.stats, wald.stats<alpha, 1)
      wald.stats=replace(wald.stats, wald.stats!=1, 0)
    }
    
    return(list(coef.mat=coef.mat, pvalues=wald.stats))
    
  }else if(response=="T3"){
    
    coef.mat=as.data.frame(summary(lm(data[,3]~., data=data[,-3]))$coefficients)
    #test.wald=coef.mat$`t value`[which(row.names(coef.mat)==colnames(data)[3])]
    wald.p.T1=coef.mat$`Pr(>|t|)`[which(row.names(coef.mat)==colnames(data)[1])]
    wald.p.T2=coef.mat$`Pr(>|t|)`[which(row.names(coef.mat)==colnames(data)[2])]
    
    wald.stats=c(wald.p.T1, wald.p.T2)
    names(wald.stats)=c("wald.p.T1","wald.p.T2")
    #returns indicator values for significant tests
    if(return.indicator==TRUE){
      wald.stats=replace(wald.stats, wald.stats<alpha, 1)
      wald.stats=replace(wald.stats, wald.stats!=1, 0)
    }
    
    return(list(coef.mat=coef.mat, pvalues=wald.stats))
  }
  
}

#----------------------Wrapper-Function-for-standard-regression--------------------
#for the alternative case when the minor variant is not "rare" and for all
#trios not including the variant
reg.without.variant=function(data=NULL, return.indicator = TRUE, alpha = 0.05){
  
  A=as.data.frame(matrix(0, nrow = 3, ncol = 3))
  row.names(A)=colnames(data[,1:3])
  colnames(A)=colnames(data[,1:3])
  
  resp.T1=do.reg.wo.var(data=data, response="T1", return.indicator = TRUE, alpha = 0.05)
  resp.T2=do.reg.wo.var(data=data, response="T2", return.indicator = TRUE, alpha = 0.05)
  resp.T3=do.reg.wo.var(data=data, response="T3", return.indicator = TRUE, alpha = 0.05)
  
  A[1,2:3]=resp.T1$pvalues
  A[2,1]=resp.T2$pvalues[1]
  A[2,3]=resp.T2$pvalues[2]
  A[3,1:2]=resp.T3$pvalues
  
  return(A)
  
}


#------------------------------Classifying-Function-----------------------------------
#a function which classifies each trio with a variant

class.Vtrios=function(pvalues=NULL){
  
  #define the model cases by regression results:
  
  M0=rbind.data.frame(c(0,0,0,1,0), c(0,1,0,0,1))
  M1=c(1,1,1,0,1)
  M2=c(1,1,1,1,0)
  M3=c(0,1,0,1,0)
  M4=c(1,1,1,1,1)
  
  cases=rbind.data.frame(M0,M1,M2,M3,M4)
  row.names(cases)=c("M0.1","M0.2","M1","M2","M3","M4")
  colnames(cases)=c("b11","b12","b21", "b22", "Cor.V.T2")
  
  
  if(all.equal(pvalues, M0[1,], check.attributes=FALSE)){
    
    
    
  }else if(all.equal(pvalues, M0[2,], check.attributes=FALSE)){
    
    
  }else if(all.equal(pvalues, M1, check.attributes=FALSE)){
    
    
  }else if(all.equal(pvalues, M2, check.attributes=FALSE)){
    
    
  }else if(all.equal(pvalues, M3, check.attributes=FALSE)){
    
    
  }else if(all.equal(pvalues, M4, check.attributes=FALSE)){
    
    
  }else{
    
    return("No matching model types")
    
  }
  
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
                                            V=V,
                                            U=U)
  
  trio.list.wo.variant=allocate.all.trios(data=data, 
                                          combs=combs.Tnodes, 
                                          with.variant=FALSE, 
                                          V=V,
                                          U=U)
  
  #get the graph skeleton
  if(is.null(U)){
    A=get.skel(data=cbind.data.frame(V,data), p=p)
  }else{
    A=get.skel(data=cbind.data.frame(V,data,U), p=p)
  }
  
  #calculate minor variant frequency
  variant.freq=get.freq(V=V, type=variant.type)
  
  if(variant.freq<gamma){
    
    out.coeffs=lapply(trio.list.with.variant, reg.with.variant, permuted=TRUE, nperms=m, return.indicator=TRUE)
    out.coeffs2=lapply(trio.list.wo.variant, reg.without.variant, return.indicator = TRUE, alpha = 0.05)
    
  }else{
    
    out.coeffs=lapply(trio.list.with.variant, reg.with.variant, permuted=FALSE, nperms=m, return.indicator=TRUE)
    out.coeffs2=lapply(trio.list.wo.variant, reg.without.variant, return.indicator = TRUE, alpha = 0.05)
    
  }
  
  
}



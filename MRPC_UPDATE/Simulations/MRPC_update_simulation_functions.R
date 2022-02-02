#These are the funtions to preform the trio-specific analysis from Section 1.1-1.2 of the pseudocode


####################################################################
#function To Get Minor Variant Frequency
get.freq=function(V=NULL){
  #used in step 2
  #remove missing values
  V=na.omit(V)
  
  #count the alternative alleles
  alternative = sum(V)/(2*length(V))
  reference=1-alternative
  
  #return the minor allele freq.
  return(min(alternative, reference))
  
}



####################################################################
#A function which preforms the standard regressions from step 1
Reg=function(data=NULL, verbose=FALSE){
  
  #data should be a n x (3+g) matrix with the variant in the first column
  #preallocate p and t vectors
  pvals=NULL
  tvals=NULL
  
  for(i in 2:3){
    
    #preform the regressions in step 1
    model=lm(data[,i]~., data = data[,-i])
    if(verbose==TRUE){print(summary(model))}
    coefs=as.data.frame(summary(model)$coefficients)
    pvals=append(pvals, coefs$`Pr(>|t|)`[-1])
    tvals=append(tvals, coefs$`t value`[-1])
    
  }
  
  #name them according to coefficient index
  names(pvals)=c("p11","p21","p12","p22")
  names(tvals)=c("tobs11","tobs21","tobs12","tobs22")
  
  #output pvalues and t-stats for b11,b21,b12,b22
  return(pt.list=list(pvals=pvals,tvals=tvals))
  
}



####################################################################
#-------------------Helper-function-for-permuted-reg---------------
PermReg.help.fn=function(perm.map=NULL, V=NULL, T1=NULL, T2=NULL, coln=NULL, response=NULL){
  
  #written for parallelization
  #permute
  if(response=="T2"){
    new.data=cbind.data.frame(V, T1[perm.map], T2)
  }else{
    new.data=cbind.data.frame(V, T1, T2[perm.map])
  }
  
  if(isFALSE(is.null(coln))){colnames(new.data)=coln}
  
  #run regression
  if(response=="T2"){
    
    #section 1.2 step 1.1 eqn (4)
    coef.mat=as.data.frame(summary(lm(new.data[,3]~., data=new.data[,-3]))$coefficients)
    wald.stat=coef.mat$`t value`[which(row.names(coef.mat)==colnames(new.data)[2])]
    #print(coef.mat)
    
  }else{
    
    #section 1.2 step 1.1 eqn (3)
    coef.mat=as.data.frame(summary(lm(new.data[,2]~., data=new.data[,-2]))$coefficients)
    wald.stat=coef.mat$`t value`[which(row.names(coef.mat)==colnames(new.data)[3])]
    #print(coef.mat)
    
  }
  
  
  #return the wald stat: section 1.2 step 1.2
  return(wald.stat)
  
}

####################################################################
#this function uses PermReg.helper.fn() to preform the permuted regression
#analysis for rare variants
PermReg=function(trio=NULL, t.obs21=NULL, t.obs22=NULL, p11=NULL, p12=NULL, m=NULL){
  
  #preallocate a matrix of indicies ranging from 1:sample size
  #we will shuffle these numbers later to get the permutations within genotype
  #need 2 for each regression: Section 1.2 steps 1-1.2 and eqn (3) and (4)
  mediator_perm1=matrix(c(1:dim(trio)[1]), nrow=dim(trio)[1], ncol = m)
  mediator_perm2=matrix(c(1:dim(trio)[1]), nrow=dim(trio)[1], ncol = m)
  
  #preallocate all permutations
  #shuffle within each genotype
  for (j in 0:2) {
    ind <- which(trio[,1] == j)
    if (length(ind) > 1) {
      mediator_perm1[ind, ] <- apply(mediator_perm1[ind, ], 2, sample)
      mediator_perm2[ind, ] <- apply(mediator_perm2[ind, ], 2, sample)
    }
  }
  
  #preforms all permutations of eqn (3) in parallel
  #outputs Theta21
  Theta21=apply(mediator_perm1, 2, PermReg.help.fn, 
                V=trio[,1], 
                T1=trio[,2], 
                T2=trio[,3], 
                coln=colnames(trio), 
                response="T1")
  #outputs Theta22
  Theta22=apply(mediator_perm2, 2, PermReg.help.fn, 
                V=trio[,1], 
                T1=trio[,2], 
                T2=trio[,3], 
                coln=colnames(trio), 
                response="T2")
  
  #Step 2.1 - calculating the nominal p-values using Theta21 and Theta22
  
  nominal.p21=2 * (1 - pnorm(abs((t.obs21 - mean(Theta21))/sd(Theta21))))
  nominal.p22=2 * (1 - pnorm(abs((t.obs22 - mean(Theta22))/sd(Theta22))))
  
  #concat
  pvals=c(p11, nominal.p21, p12, nominal.p22)
  #return pvalue vector
  return(pvals)
}




####################################################################
#a wrapper function for get.freq(), Reg(), and PermReg() to infer the trio
#combines the functions from sections 1.1-1.2
infer.trio=function(trio=NULL, gamma=0.05, alpha=0.01, nperms=1000, verbose=FALSE){
  
  #preallocate indicator vectors
  xp=NULL
  rp=NULL
  
  #preform the standard regressions and outputs t-stat and p-values
  #input is a trio with the variant in the first column
  #step 1
  pt.out=Reg(data = trio, verbose=verbose)
  

  
  #check the frequency of the minor allele using get.freq()
  #step 2
  minor=get.freq(V=trio[,1])
  
  if(minor<gamma){
    #preform permuted regression for rare variants
    pvals=PermReg(trio = trio, 
                  t.obs21 = pt.out$tvals[2], 
                  t.obs22 = pt.out$tvals[4], 
                  p11 = pt.out$pvals[1],
                  p12 = pt.out$pvals[3],
                  m = nperms)
    
  }else{
    #else return the pvals from standard reg.
    pvals=pt.out$pvals
    
  }
  
  #---steps 3-4
  #convert xp to indicator vector
  #section 1.1 step
  xp=ifelse(pvals<alpha, 1, 0)
  
  #preform marginal tests
  cors=c(cor.test(trio[,1], trio[,3])$p.value, cor.test(trio[,1], trio[,2])$p.value)
  #convert to indicator vector
  rp=ifelse(cors<alpha, 1, 0)
  #combine all useful stats - add indicator 
  all.stats=c(append(xp, rp), append(pvals, cors), minor)
  
  names(all.stats)=c("b11","b21", "b12","b22", "V1:T2", "V1:T1", "pb11", 
                     "pb21", "pb12","pb22","pV1:T2","pV1:T1", "Minor.freq")
  
  return(all.stats)
}


####################################################################
#A function to classify each indicator vector returned by inf.trio()
#section 1.1 steps 4-5
class.vec=function(vec=NULL){


  M0=matrix(c(1,0,0,0,0,0,1,0), nrow = 2, ncol = 4, byrow = T)
  M1=matrix(c(1,1,0,1,0,1,1,1), nrow = 2, ncol = 4, byrow = T)
  M2.M4=c(1,1,1,1)
  M3=c(1,0,1,0)

  ind.mat=rbind.data.frame(M0,M1,M3,M2.M4)
  row.names(ind.mat)=c("M0.1","M0.2","M1.1","M1.2","M3","M2/M4")

  which.mod=row.names(ind.mat)[row.match(vec[1:4], ind.mat)]
  #print(which.mod)
  if(is.na(row.match(vec[1:4], ind.mat))){ct="Other"}
  
  else if(which.mod=="M2/M4"){
    
    rp=vec[5:6]
    if(sum(rp)==2){ct="M4"}
    
    else if(sum(rp==c(0,1))==2){ct="M2.1"}
    
    else if(sum(rp==c(1,0))==2){ct="M2.2"}
    
    else{ct="Other"}
    
  }else{ct=which.mod}

  return(ct)

}


####################################################################\
#Version 2
#A function to classify each indicator vector returned by inf.trio()
#section 1.1 steps 4-5
# class.vec=function(vec=NULL){
#   
#   number.of.pos=sum(vec[1:4])
#   M0=matrix(c(1,0,0,0,0,0,1,0), nrow = 2, ncol = 4, byrow = T)
#   #row.names(M0)=c("M0.1","M0.2")
#   M1=matrix(c(1,1,0,1,0,1,1,1), nrow = 2, ncol = 4, byrow = T)
#   #row.names(M1)=c("M1.1","M1.2")
#   M2.M4=c(1,1,1,1)
#   M3=c(1,0,1,0)
#   
#   ind.mat=rbind.data.frame(M0,M1,M2.M4,M3)
#   row.names(ind.mat)=c("M0.1","M0.2","M1.1","M1.2","M2/M4","M3")
#   which.mod=row.names(ind.mat)[row.match(vec[1:4], ind.mat)]
#   
#   if(which.mod=="M3"){ct="M3"}
#   
#   else if(which.mod=="M0.1" | which.mod=="M0.2" & sum(vec[5:6])==1){
#     
#     ct=which.mod
#     
#   }
#   
#   else if(which.mod=="M1.1" | which.mod=="M1.2" & sum(vec[5:6])==2){
#     ct=which.mod
#   }
#   
#   else if(which.mod=="M2/M4" & sum(vec[5:6])==2){ct="M4"}
#   
#   else if(which.mod=="M2/M4" & sum(vec[5:6])==1){
#     
#     ct=ifelse(vec[5]==0, "M2.1", "M2.2")
#     
#   }else{ct="Other"}
#     
# 
#   
#   return(ct)
#   
# }


####################################################################\
#Version 3
#A function to classify each indicator vector returned by inf.trio()
#section 1.1 steps 4-5
# class.vec=function(vec=NULL){
#   
#   number.of.pos=sum(vec[1:4])
#   M0=matrix(c(1,0,0,0,0,0,1,0), nrow = 2, ncol = 4, byrow = T)
#   #row.names(M0)=c("M0.1","M0.2")
#   M1=matrix(c(1,1,0,1,0,1,1,1), nrow = 2, ncol = 4, byrow = T)
#   #row.names(M1)=c("M1.1","M1.2")
#   M2.M4=c(1,1,1,1)
#   M3=c(1,0,1,0)
#   
#   ind.mat=rbind.data.frame(M0,M1,M2.M4,M3)
#   row.names(ind.mat)=c("M0.1","M0.2","M1.1","M1.2","M2/M4","M3")
#   which.mod=row.names(ind.mat)[row.match(vec[1:4], ind.mat)]
#   
#   if(which.mod=="M3"){ct="M3"}
#   
#   else if(which.mod=="M0.1" | which.mod=="M0.2" & sum(vec[5:6])==1){
#     
#     ct=which.mod
#     
#   }
#   
#   else if(which.mod=="M1.1" | which.mod=="M1.2" & sum(vec[5:6])==2){
#     ct=which.mod
#   }
#   
#   else if(which.mod=="M2/M4" & sum(vec[5:6])==2){ct="M4"}
#   
#   else if(which.mod=="M2/M4" & sum(vec[5:6])==1){
#     
#     ct=ifelse(vec[5]==0, "M2.1", "M2.2")
#     
#   }else{ct="Other"}
#   
#   
#   
#   return(ct)
#   
# }
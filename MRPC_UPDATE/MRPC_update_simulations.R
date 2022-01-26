
# a script to simulate multiple types of trios and their indicator vector under the 
#MRPC updated framework to determine if the current set of cases is enough to identify
#each of the model types M0,M1,M2,M3,M4

#load necessary packages
#relevant packages
library(MRPC)
library(psych)
library(prodlim)
library('MRPC', lib="/mnt/ceph/jarredk/Rpackages")
library('qvalue', lib="/mnt/ceph/jarredk/Rpackages")
library('psych', lib="/mnt/ceph/jarredk/Rpackages")
library('prodlim', lib="/mnt/ceph/jarredk/Rpackages")

# source("/mnt/ceph/jarredk/MRPC_UPDATE/MRPCgeneral.R")
# source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")


#set initial parameters for simulation
#sample sizes
n=c(50, 100, 500, 1000)
#noise in the data/SD of errors
noise=c(0.2, 0.5, 0.8, 1.2)
#frequency of the minor allele
minor.allele=c(0.1, 0.2, 0.3)
#signal strength of edge
b1.1=c(0.2, 0.4, 0.6)
b1.2=c(0.2, 0.4, 0.6)

#expand to all get all combos of simulation conditions
model.params=expand.grid(n, noise, minor.allele, b1.1, b1.2)
colnames(model.params)=c("sample.size","SD", "minor.freq", "b1.1","b1.2")
#the number of simulations for each set of conditions
sims=100

####################################################################
#a function to preform the standard regression for minor allele frequencies above gamma
regress=function(trio=NULL, alpha=0.01, verbose=FALSE){
  
  indicator=NULL
  pvals=NULL
  
  for(i in 1:2){
    
    model=lm(trio[,i+1]~., data = trio[,-(i+1)])
    if(verbose==TRUE){print(summary(model))}
    coefs=as.data.frame(summary(model)$coefficients)
    pvals=append(pvals, coefs$`Pr(>|t|)`[-1])
    
  }
  
  indicator=ifelse(pvals<alpha, 1, 0)
  
  #cors=corr.test(trio)$p[1,2:3]
  
  #indicator=append(indicator, ifelse(cors<alpha, 1, 0))
  
  #names(indicator)=c("b11","b21", "b12","b22", "V1:T1", "V1:T2")
  names(indicator)=c("b11","b21", "b12","b22")
  return(indicator)
}


####################################################################
#A function to classify each indicator vector returned by regress()
class.vec=function(vec=NULL, trio=NULL, alpha=0.01){
  
  M0=matrix(c(1,0,0,0,0,0,1,0), nrow = 2, ncol = 4, byrow = T)
  M1=matrix(c(1,1,0,1,0,1,1,1), nrow = 2, ncol = 4, byrow = T)
  M2.M4=matrix(c(1,1,1,1), nrow = 1, ncol = 4, byrow = T)
  M3=matrix(c(1,0,1,0), nrow = 1, ncol = 4, byrow = T)
  
  ind.mat=rbind.data.frame(M0,M1,M2.M4,M3)
  row.names(ind.mat)=c("M0.1","M0.2","M1.1","M1.2","M2/M4","M3")
  
  which.mod=row.names(ind.mat)[row.match(vec, ind.mat)]
  #print(which.mod)
  
  #a series of {if} statements to parse each trio into one of the 
  #5 topologies (or into the "Other" class if no topology matches)
  if(is.na(which.mod)==TRUE){ct="Other"}else{
    
    if(which.mod=="M0.1"){
      ct=ifelse(corr.test(trio)$p[1,3]>alpha, "M0", "Other")
    }else if(which.mod=="M0.2"){
      ct=ifelse(corr.test(trio)$p[1,2]>alpha, "M0", "Other")
    }else if(which.mod=="M1.1"){
      ct=ifelse(corr.test(trio)$p[1,3]<alpha, "M1", "Other")
    }else if(which.mod=="M1.2"){
      ct=ifelse(corr.test(trio)$p[1,2]<alpha, "M1", "Other")
    }else if(which.mod=="M2/M4"){
      test=ifelse(corr.test(trio)$p[1,2:3]<alpha, 1, 0)
      if(sum(test)==1){
        ct="M2"
      }else if(sum(test)==2){
        ct="M4"
      }else{
        ct="Other"
      }
    }else if(which.mod=="M3"){
      ct="M3"
    }
    
  }
  
  
  
  return(ct)
  
}

####################################################################

##----------M0-----------

store.vectors=list()
#data.sets=list()
accuracy=NULL

#looping for all combos of simulation conditions
for(i in 1:dim(model.params)[1]){
  
  #pre-allocate for inner loop
  vec.mat=as.data.frame(matrix(NA, nrow = sims, ncol=4))
  colnames(vec.mat)=c("b11", "b21", "b12", "b22")
  classes=NULL
  #all.sim=list()
  
  #run simulations
  for(j in 1:sims){
    #simulate trio under parameters
    X=SimulateData(N=model.params$sample.size[i], 
                   p=model.params$minor.freq[i], 
                   model="model0", 
                   b0.1=0, 
                   b1.1=model.params$b1.1[i], 
                   sd.1=model.params$SD[i])
    
    #store regression test indicator vector
    vec.mat[j,]=regress(X)
    #print(vec.mat[j,])
    #Store inferred class
    classes[j]=class.vec(vec = vec.mat[j,], trio = X)
    #all.sim[[j]]=X
  }
  #print(vec.mat)
  #store inner loop objects into outer loop lists
  store.vectors[[i]]=cbind.data.frame(vec.mat, inf.class=classes)
  #calculate a basic prediction accuracy for each set of simulation conditions
  accuracy[i]=length(which(store.vectors[[i]]$inf.class="M0"))/dim(store.vectors[[i]])[1]
  
  print(paste0("finished sim ", i, " of ", dim(model.params)[1]))
  
}





##----------M1-----------





##----------M2-----------



##----------M3-----------





##----------M4-----------
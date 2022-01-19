
# a script to simulate multiple types of trios and their indicator vector under the 
#MRPC updated framework to determine if the current set of cases is enough to identify
#each of the model types M0,M1,M2,M3,M4

#load necessary packages

library('MRPC', lib="/mnt/ceph/jarredk/Rpackages")
library('qvalue', lib="/mnt/ceph/jarredk/Rpackages")

source("/mnt/ceph/jarredk/MRPC_UPDATE/MRPCgeneral.R")
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
l1=cross.analyze(tissues="AdiposeSubcutaneous", save=FALSE)

trio.list=list()
for(i in 1:25){
  
  list.data=cross.regress(tissue="AdiposeSubcutaneous", trio.ind=l1$final.tables[[1]]$Trio.Num[i], mod.type="both", addis.pcs=NULL, verbose=F)
  data=list.data$GMAC[, -((dim(list.data$GMAC)[2]-2):dim(list.data$GMAC)[2])]
  data=data[,c(3,2,1, 4:dim(data)[2])]
  trio.list[[i]]=as.data.frame(data)
}

m=1000
out.coeffs=lapply(trio.list, reg.with.variant, permuted=TRUE, nperms=m, return.indicator=TRUE, alpha=0.01)


loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#load in necesarry files from workspace

load(file="/mnt/ceph/jarredk/Methyl/Wksp1.Rdata")

genos.mat=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/fakegenos.Rdata")

EM.triolist=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/completedTrios1.Rdata")

triobuildlist=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/triobuildlist.Rdata")



#function to compute correlations between SNPS 


calc.cors=function(triolist=NULL, genomat=NULL){
  
  size1=NULL
  for(i in 1:length(triobuildlist)){size1[i]=length(triobuildlist[[i]])}
  
  print(paste("number of trios to compute", paste0(sum(vec1),"X",dim(genos.mat)[2]), sep = ":"))
  
  
  
}












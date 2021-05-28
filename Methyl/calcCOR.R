
#================================================================================
#--------------------build-pieces-needed-to-run-correlations---------------------
#================================================================================


loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#load in necesarry files from workspace

load(file="/mnt/ceph/jarredk/Methyl/Wksp1.Rdata")
genos.mat=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/fakegenosBIG.Rdata")
EM.triolist=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/completedTrios1.Rdata")
triobuildlist=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/triobuildlist.Rdata")
#meta_G=read.delim(file="/mnt/ceph/megheib/M_G_data/GPL10558_GE_R2.txt", sep="")
meta_M=read.csv(file = "/mnt/ceph/megheib/M_G_data/GPL13534_M.csv")

imp.meta_M=meta_M[,c(1,12,13)]

#simulate the chromosome and position of each SNP

sim.chr=sample(c(1:22,"X","Y"), dim(genos.mat)[2], replace = TRUE)

simu.pos=function(genosmat=NULL, chr.vec=NULL){
  
  #pos=matrix(0, nrow = dim(genosmat)[2], ncol = 2)
  pos=NULL
  chr.sizes=read.csv(file = "/mnt/ceph/jarredk/Methyl/chrsizes.csv")
  
  for(i in 1:length(chr.vec)){
    
    idx=match(chr.vec[i],chr.sizes[,1])
    start.pos=sample(c(1000000:(chr.sizes[idx,2]-1000000)), 1)
    #end.pos=start.pos+1000000
    #print(start.pos)
    #print(end.pos)
    pos[i]=start.pos
    
  }
  
  return(simu.pos)
  
}

sim.chrpos=simu.pos(genosmat = genos.mat, chr.vec = sim.chr)
simulated.meta=cbind.data.frame(sim.chr, simu.pos)

#save(simulated.meta, file = "/mnt/ceph/jarredk/Methyl/fakegenoMeta.Rdata")
simulated.meta=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/fakegenoMeta.Rdata")

#acquire the "in-use" probe ID's from the completed trio list

probenames=list()
cn=function(X){ v=colnames(X)[2]; return(v)}

for(i in 1:length(EM.triolist)){
  
  if(is.null(unlist(lapply(EM.triolist[[i]], cn)))){
    
    probenames[[i]]=NA
    
  }else{
    
    probenames[[i]]=unlist(lapply(EM.triolist[[i]], cn))
    
  }
}

Mprobenames.final=na.omit(unlist(probenames))

#save(Mprobenames.final, file = "/mnt/ceph/jarredk/Methyl/Used.Mprobes.final.Rdata")
Mprobenames.final=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/Used.Mprobes.final.Rdata")


#retreive relevant metadata for the "In-Use" methylation probes:

matched.meta_M=match(Mprobenames.final, imp.meta_M[,1])
imp.meta_M.final=imp.meta_M[matched.meta_M,]




















#=========================================================================================================
#-------------------------------------calculate-the-correlations------------------------------------------
#=========================================================================================================
#function to compute correlations between SNPS based on locations within each chromosome
#load in the necessary information:

genos.mat=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/fakegenosBIG.Rdata")
EM.triolist=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/completedTrios1.Rdata")
triobuildlist=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/triobuildlist.Rdata")
simulated.meta=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/fakegenoMeta.Rdata")
Mprobenames.final=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/Used.Mprobes.final.Rdata")
meta_M=read.csv(file = "/mnt/ceph/megheib/M_G_data/GPL13534_M.csv")
#retreive relevant metadata for the "In-Use" methylation probes:
imp.meta_M=meta_M[,c(1,12,13)]
matched.meta_M=match(Mprobenames.final, imp.meta_M[,1])
imp.meta_M.final=imp.meta_M[matched.meta_M,]



calc.cors=function(EMtriolist=NULL, GMInfo=NULL, genomat=NULL, genoInfo=NULL){
  
  #size1=NULL
  #for(i in 1:length(triobuildlist)){size1[i]=length(triobuildlist[[i]])}
  #print(paste("number of trios to compute", paste0(sum(vec1),"X",dim(genos.mat)[2]), sep = ":"))
  
  
  
}












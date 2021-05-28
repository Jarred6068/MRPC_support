
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
  pos=rep(0, 4000000)
  chr.sizes=read.csv(file = "/mnt/ceph/jarredk/Methyl/chrsizes.csv")
  
  for(i in 1:length(chr.vec)){
    
    idx=match(chr.vec[i],chr.sizes[,1])
    start.pos=runif(1, 1000000, chr.sizes[idx,2]-1000000)
    #end.pos=start.pos+1000000
    #print(start.pos)
    #print(end.pos)
    pos[i]=start.pos
    
  }
  
  return(pos)
  
}

sim.chrpos=simu.pos(genosmat = genos.mat, chr.vec = sim.chr)
simulated.meta=cbind.data.frame(sim.chr, simu.chrpos)

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

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

load(file="/mnt/ceph/jarredk/Methyl/Wksp1.Rdata")
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



#==================helper_function_1===============================
# fn1=function(X,Y){ match(Y, colnames(X))}
# lapply(a, fn1, "a.1")



prep.cors=function(EMtriolist=NULL, GMInfo=NULL, genomat=NULL, genoInfo=NULL){
  
  #size1=NULL
  #for(i in 1:length(triobuildlist)){size1[i]=length(triobuildlist[[i]])}
  #print(paste("number of trios to compute", paste0(sum(vec1),"X",dim(genos.mat)[2]), sep = ":"))
  
  #modify existing triolist
  EMtriolist.mod=vector("list",length(EMtriolist))
  
  for(i in 1:length(EMtriolist)){
    
    EMtriolist.mod[[i]]=data.frame(EMtriolist[[i]])
    
  }
  
  #chrs=c(1:22,"X","Y")
  chrs=c("1","2")
  probes.by.chr=vector("list", length(chrs))
  genos.by.chr=vector("list", length(chrs))
  
  for(i in 1:length(chrs)){
    
    snps.in.range=list()
    chr.geno.matched=genoInfo[which(genoInfo[,1]==chrs[i]),]
    chr.GM.matched=GMInfo[which(GMInfo[,2]==chrs[i]),]
    
    
    
    for(j in 1:dim(chr.GM.matched)[1]){
      
      #index the snps that are close to each probe
      snps.in.range[[j]]=which(abs(chr.GM.matched[j,3]-chr.geno.matched[,2])<1000000)
      
      #print(probes.in.range[[i]])
    }
    
    #store the snps corresponding to probes on each chromosome in master list
    names(snps.in.range)=chr.GM.matched[,1]
    probes.by.chr[[i]]=snps.in.range
    #genos.by.chr[[i]]=which(genoInfo[,1]==chrs[i])
    
  }
  
  #now calculated the correlations from the indexes
  
  
  names(probes.by.chr)=chrs
  
  return(list(Mprobes=probes.by.chr, genos=genos.by.chr, EM.triolist.Mod=EMtriolist.mod))
  
}

try1=prep.cors(EMtriolist=EM.triolist, GMInfo=imp.meta_M.final, genomat=genos.mat, genoInfo=simulated.meta)










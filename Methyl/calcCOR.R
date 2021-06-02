
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

names1=paste0("SNP",c(1:4000000), sim.chr)

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
#-------------------------------prep-data----------------------------------
loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#load in previously cleaned data
load(file="/mnt/ceph/jarredk/Methyl/Wksp1.Rdata")
genos.mat=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/fakegenosBIG.Rdata")
#EM.triolist=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/completedTrios1.Rdata")
#triobuildlist=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/triobuildlist.Rdata")
simulated.meta=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/fakegenoMeta.Rdata")
Mprobenames.final=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/Used.Mprobes.final.Rdata")

#load in methylation meta-data
meta_M=read.csv(file = "/mnt/ceph/megheib/M_G_data/GPL13534_M.csv")
#load in expression meta-data from biomart
meta_E=read.table(file = "/mnt/ceph/jarredk/Methyl/mart_export.txt", sep = ",", header = T)

#retrieve relevant metadata for the "In-Use" methylation probes:
imp.meta_M=na.omit(cbind(meta_M[,c(1,15)], as.numeric(meta_M[,16])))
matched.meta_M=match(Mprobenames.final, imp.meta_M[,1])
imp.meta_M.final=imp.meta_M[matched.meta_M,]
#retrieve relevant metadata for the "In-Use" Expression probes:
matched.meta_E=match(express.genenames, meta_E[,1])
imp.meta_E.final=meta_E[na.omit(matched.meta_E), ]
#filter out columns (genes) that didn't have chr position information according to BioMart
testexpress=express.aligned[,-attr(na.omit(matched.meta_E), "na.action")]
testexpress.gn=express.genenames[-attr(na.omit(matched.meta_E), "na.action")]
#-------------------------------------------------------------------------
#=========================helper_function_1===============================

prep.cors=function(GMInfo=NULL, GEInfo=NULL, genomat=NULL, genoInfo=NULL, chrs=NULL){
  #SYNTAX:
  
  #GMInfo -- Gene Methylation matrix additional info in the form
                    # [ probe_name  Chr  start_position  ]   
  
  #GEInfo -- Gene Expression matrix of additional info in the form
                    # [ Gene_name  Chr  Start_position  ]
  
  #genomat -- the matrix/df of genotype data (with columnnames)
  
  #genoInfo -- SNP matrix/df of additional info in the form
                   # [ Chr  Start_position  ]
  
  #Chrs -- passed from calc.cors
  

  Mprobes.by.chr=vector("list", length(chrs))
  Eprobes.by.chr=vector("list", length(chrs))
  genos.by.chr=vector("list", length(chrs))
  
  for(i in 1:length(chrs)){
    
    snps.in.range=list()
    snps.in.range2=list()
    chr.geno.matched=genoInfo[which(genoInfo[,1]==chrs[i]),]
    chr.GM.matched=GMInfo[which(GMInfo[,2]==chrs[i]),]
    chr.GE.matched=GEInfo[which(GEInfo[,2]==chrs[i]), ]
    #print(chr.GE.matched)
    
    for(j in 1:dim(chr.GM.matched)[1]){
      
      #index the snps that are close to each Methylation probe
      snps.in.range[[j]]=which(abs(chr.GM.matched[j,3]-chr.geno.matched[,2])<1000000)
      #print(snps.in.range[[j]])
    }
    
    for(j in 1:dim(chr.GE.matched)[1]){
      
      #index the snps that are close to each Expression probe
      snps.in.range2[[j]]=which(abs(chr.GE.matched[j,3]-chr.geno.matched[,2])<1000000)
      #print(snps.in.range2[[j]])
    }
    
    #store the snps corresponding to probes on each chromosome in master list
    names(snps.in.range)=chr.GM.matched[,1]
    names(snps.in.range2)=chr.GE.matched[,1]
    Mprobes.by.chr[[i]]=snps.in.range
    Eprobes.by.chr[[i]]=snps.in.range2
    #genos.by.chr[[i]]=which(genoInfo[,1]==chrs[i])
    
  }
  
  #now calculated the correlations from the indexes
  
  names(Mprobes.by.chr)=paste0("chr",chrs)
  names(Eprobes.by.chr)=paste0("chr",chrs)
  return(list(Mprobes=Mprobes.by.chr, Eprobes=Eprobes.by.chr))
  
}

#test code
try1=prep.cors(GMInfo=imp.meta_M.final, GEInfo=imp.meta_E.final, genomat=genos.mat, genoInfo=simulated.meta, chrs = c("1"))


#==========================================function2==============================================

calc.cors=function(mmat=NULL, emat=NULL, gmat=NULL, GMInfo=NULL, GEInfo=NULL, genoInfo=NULL, ge.names=NULL, chrs=NULL){
  #SYNTAX:
  
  #mmat - the methylation matrix/df with column-names as probe ID's
  
  #emat - the expression matrix (aligned with BIOmart meta-data in GEInfo)
  
  #GMInfo -- Gene Methylation matrix additional info in the form
  # [ probe_name  Chr  start_position  ]   
  
  #GEInfo -- Gene Expression matrix of additional info in the form
  # [ Gene_name  Chr  Start_position  ]
  
  #genomat -- the matrix/df of genotype data (with column-names as SNP ID's)
  
  #genoInfo -- SNP matrix/df of additional info in the form
  # [ Chr  Start_position  ]
  
  #ge.names -- the gene names/ID's corresponding to the columns of 'emat' 
  
  #Chrs -- The chromosome(s) for which data is to be searched and saved (value between "1:22", "X", "Y")
  
  
  #run prep cors to get indexes of SNPs corresponding to each M and E probe
  prepped1=prep.cors(GMInfo=GMInfo, GEInfo = GEInfo, genomat=gmat, genoInfo=genoInfo, chrs=chrs)
  
  #pre-allocation and naming of storage lists
  methyl.cors.list=vector("list", length(prepped1$Mprobes))
  express.cors.list=vector("list", length(prepped1$Eprobes))
  names(methyl.cors.list)=names(prepped1$Mprobes)
  names(express.cors.list)=names(prepped1$Eprobes)
  
  #calculate correlations with methyl probes
  
  for(k in 1:length(chrs)){
    
    gmat.new=gmat[,which(genoInfo[,1]==chrs[k])]
    
    for(i in 1:length(prepped1$Mprobes)){
      
      snp.cors=list()
      cols=match(names(prepped1$Mprobes[[i]]), colnames(mmat))
      
      methyl.mat=mmat[,cols]
      
      for(j in 1:dim(methyl.mat)[2]){
        
        cor_j=cor(methyl.mat[,j], gmat.new[, prepped1$Mprobes[[i]][[j]] ], use = "complete.obs")
        names(cor_j)=colnames(gmat.new[, prepped1$Mprobes[[i]][[j]] ])
        snp.cors[[j]]=cor_j
      }
      
      names(snp.cors)=names(prepped1$Mprobes[[i]])
      
      methyl.cors.list[[i]]=snp.cors
    }
    
    #calculate correlations with Expression probes
    
    for(i in 1:length(prepped1$Eprobes)){
      
      snp.cors2=list()
      cols2=match(names(prepped1$Eprobes[[i]]), ge.names)
      express.mat=emat[,cols2]
      
      for(j in 1:dim(express.mat)[2]){
        
        cor_j2=cor(express.mat[,j], gmat.new[, prepped1$Eprobes[[i]][[j]] ], use = "complete.obs")
        names(cor_j2)=colnames(gmat.new[, prepped1$Eprobes[[i]][[j]] ])
        snp.cors2[[j]]=cor_j2
      }
      
      names(snp.cors2)=names(prepped1$Eprobes[[i]])
      
      express.cors.list[[i]]=snp.cors2
    }
    
    
  }
  
  
  
  
  #calculate correlations with methyl probes
  return(list(Mcors=methyl.cors.list, Ecors=express.cors.list))
  
}

#========================================END_Funtion===========================================

#test code

trycors1=calc.cors(mmat = methyl.resids2,
                   emat = testexpress,
                   gmat = genos.mat, 
                   GMInfo = imp.meta_M.final, 
                   GEInfo = imp.meta_E.final, 
                   genoInfo = simulated.meta, 
                   ge.names = testexpress.gn,
                   chrs=c("1"))




#calculate for all available/usable probes and save by chromosome (test example)

cors1=calc.cors(mmat = methyl.resids2,
                   emat = testexpress,
                   gmat = genos.mat, 
                   GMInfo = imp.meta_M.final, 
                   GEInfo = imp.meta_E.final, 
                   genoInfo = simulated.meta, 
                   ge.names = testexpress.gn,
                   chrs=c("1","2"))

save(cors1, file = "/mnt/ceph/jarredk/Methyl/EXcorLists/cors_chr1_chr2.Rdata")
cors1=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/EXcorLists/cors_chr1_chr2.Rdata")








































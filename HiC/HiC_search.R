
#read in ADDIS_Post_Analysis_processing files
source("/mnt/ceph/jarredk/ADDIS_verify/ADDIS_Post_Analysis_processing.R")
source("/mnt/ceph/jarredk/ADDIS_verify/ADDIS_rerun.R")

#A simple helper function for easier reading in of .Rdata files with a specific name
#================================================Helper_Function============================================

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#==============================================Helper_Function==============================================

find_trans=function(TransGeneName=NULL){
  
  tissue.names=read.csv("/mnt/ceph/jarredk/AddisReRunFiles/tissuenames.csv", header = T)

    meta2=read.delim(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissue.names[1,2],
                    "_AllPC/",tissue.names[1,2],".genepos.V8.txt",sep="")) 
      
      TGattr=subset(meta2, gene_id==TransGeneName)
  
  return(TGattr)
  
}

#==================================Begin_Helper_Function_2==================================================
get_trio_attr=function(trio.index=NULL, tissue.name="CellsEBVtransformedlymphocytes", verbose=TRUE){
  #SYNTAX: 
  #trio.index -- the index for the specific cis/trans mediated trio (M1's)
  #index.tissue -- currently only one value ='LCL'
  #a simple checking mechanism
  
  tissue.names=read.csv("/mnt/ceph/jarredk/AddisReRunFiles/tissuenames.csv", header = T)
  name=subset(tissue.names, tissue.name1==tissue.name)
  
  # load necessary files for given tissue
    
  triodata = loadRData(fileName = paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/", 
                                          tissue.name,"_AllPC/data.snp.cis.trans.final.",tissue.name,".V8.unique.snps.RData", sep = ""))
  
  meta=loadRData(fileName = paste("/mnt/ceph/jarredk/AddisReRunFiles/", name[3], ".eGenes.V8.unique.RData", sep = ""))
  
  meta2=read.delim(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissue.names[1,2],
                         "_AllPC/",tissue.names[1,2],".genepos.V8.txt",sep="")) 
    
  #space allocation
  attri.mat.cis=as.data.frame(matrix(0, nrow=length(trio.index), ncol=dim(meta)[2]))
  attri.mat.trans=as.data.frame(matrix(0, nrow=length(trio.index), ncol=4))
  snps=NULL
  #retrieve attributes for all genes in all given trios: 
  for(i in 1:length(trio.index)){
      
    #find trio and extract
    gene.name=NULL
    trio=triodata[,(trio.index[i]-1)*3+(1:3)]
    gene.name=colnames(trio)[-1]
    snps[i]=colnames(trio)[1]
    
    if(verbose==TRUE){
      print(paste(i,") finding attributes for ", trio.index[i],"...", sep = ""))
    }
      
    #get cis-trans gene attributes
    TGattr=find_trans(TransGeneName=gene.name[2])
    attri.mat.cis[i,] = subset(meta, gene_id == gene.name[1])
    attri.mat.trans[i,] = TGattr
      
  }
  
  colnames(attri.mat.trans)=colnames(meta2)
  colnames(attri.mat.cis)=colnames(meta)
  attri.mat.cis$chr = gsub("chr","", attri.mat.cis$chr)
  attri.mat.trans$trio.idx=trio.index
  attri.mat.cis$trio.idx=trio.index
  data.list=list(trans=attri.mat.trans, cis=attri.mat.cis)
  
  return(list(Attributes=data.list, SNPS=snps))
  
  
}

#==================================helper_Function_3=======================================================

extract_hic=function(fileName=NULL, chrs=c("1","1"), resol=10000, package.path="/mnt/ceph/jarredk/Rpackages", savefile=FALSE){
  #SYNTAX:
  #fileName -- the path name to the .hic file 
  #chrs -- (option1) the chromosomes for which chromatin interactions are to be measured: if both numbers are the same chromosome 
  #the sparse upper matrix from the .hic returned is interchromosomal interactions: (option2) reads from specific regions can be 
  #extracted: E.G.) -- Extract reads between 1MB and 7.5MB on chromosome 1 using c("1:1000000:7400000", "1:1000000:7400000")
  #package.path -- the path to your R-package repo
  #savefile -- logical indicating if the extracted file should be saved
  
  library(strawr, lib=package.path)
  
  HiC.df=straw("observed", "NONE", fileName, chrs[1], chrs[2], "BP", resol)
  
  if(savefile==TRUE){
    write.table(HiC.df, file = paste(fileName, ".txt", sep=""), sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  
  return(HiC.df)
  
}


#=============================================Resampling_Function====================================================

Resample_interactions=function(filePath=NULL, chrs=c("1","1"), res=10000, search.size=100000, resamples=10000, 
                               verbose=FALSE, plot.hist=FALSE){
  
  reads=rep(0, resamples)
  sampler=NULL
  CI=NULL
  i=1
  data.hic2=extract_hic(fileName = filePath, 
                        chrs = chrs,
                        resol=res)
  
  bounds=c(min(data.hic2$x), max(data.hic2$x), min(data.hic2$y), max(data.hic2$y))
  
  if(verbose==TRUE){
    print(bounds)
  }

  for( i in 1:length(reads)){
    
    gene1=runif(2, min = bounds[1]+search.size, max = bounds[2]-search.size)
    gene2=runif(2, min = bounds[3]+search.size, max = bounds[4]-search.size)
    
    #establish snp - gene search "box"
    left.pos=paste(chrs[1],":",
                   gene1[1]-search.size,":",
                   gene1[1]+search.size, sep = "")
  
    right.pos=paste(chrs[2],":",
                    gene2[2]-search.size,":",
                    gene2[2]+search.size, sep = "")
  
  
    hic.obj=extract_hic(fileName=filePath, chrs = c(left.pos, right.pos), resol = res)
    
    #count interactions
    reads[i]=ifelse(empty(as.data.frame(hic.obj))==TRUE, NA, sum(hic.obj$counts))
    
    # if(is.na(sampler[i])!=TRUE){
    #   reads[i]=sampler[i]
    # }
    # 
    # i=i+1
    
  }
  
  
  if(plot.hist==TRUE){
    
    H1=hist(reads)
    plot(H1)
    
  }
  
  reads=na.omit(reads)
  
  X.bar=mean(reads)
  X.sd=sd(reads)
  lower.limit=X.bar-abs(qnorm(0.025))*X.sd
  upper.limit=X.bar+abs(qnorm(0.025))*X.sd
  CI=c(lower.limit, X.bar, upper.limit)
  names(CI)=c("lower.limit", "mean", "upper.limit")
  
  
  if(plot.hist==TRUE){
    return(list(resampled.totals=reads,rawdata=sampler, confInterval=CI, Hist.obj=H1))
  }else{
    return(list(resampled.totals=reads, confInterval=CI))
  }
  
  
}


#==================================Check_Interactions_Function=======================================================

# a simple function which incorporates the HiC data in .hic format to check if chromatin interactions exist for trios
# classed as trans or cis mediated models (M1's) obtained from post-processing analysis. 

interaction_check=function(hic.filename=NULL, trios=NULL, resolution=10000, search.size=100000, tiss="CellsEBVtransformedlymphocytes",
                           verbose=TRUE){
  #SYNTAX:
  #hic.filename -- the path to the desired hic.file
  #trio.attr -- the attributes information from "get_trio_attr()" 
  #resolution -- the BP resolution desired and to be passed to "extract_hic()"
  
  
  trio.attr=get_trio_attr(trio.index=trios, tissue.name = tiss)
  
  cis.data=trio.attr$Attributes$cis
  trans.data=trio.attr$Attributes$trans
  
  idx = length(trios)
  reads = NULL
  averages = NULL
  p.values = NULL
  
  hic.extent=as.data.frame(matrix(0, nrow = idx, ncol = 4))
  colnames(hic.extent)=c("lower_x", "upper_x", "lower_y", "upper_y")
  
  for(i in 1:idx){
    #establish snp - gene search "box"
    left.pos=paste(cis.data$chr[i],":",
                   cis.data$variant_pos[i]-search.size,":",
                   cis.data$variant_pos[i]+search.size, sep = "")
    
    #establish trans gene search "box"
    right.pos=paste(trans.data$chr[i],":",
                    trans.data$left[i]-search.size,":",
                    trans.data$right[i]+search.size, sep = "")
    
    
    if(verbose==TRUE){
      print(left.pos)
      print(right.pos)
    }
    
    
    #save search upper and lower bounds
    hic.extent[i,]=c(cis.data$variant_pos[i]-search.size,
                     cis.data$variant_pos[i]+search.size,
                     trans.data$left[i]-search.size,
                     trans.data$right[i]+search.size)
    
    #extract data from .hic based on search parameters
    data.hic=extract_hic(fileName = hic.filename, 
                         chrs = c(left.pos, right.pos),
                         resol=resolution)
    
    #count interactions in observed
    reads[i]=ifelse(empty(as.data.frame(data.hic))==TRUE, NA, sum(data.hic$counts))
    
    
    #resample the data to obtain a probability distribution by MC integration
    if(is.na(reads[i])==TRUE){
      averages[i]=NA
      p.values[i]=NA
    }else{

      RS=Resample_interactions(filePath = hic.filename,
                              chrs=c(paste(cis.data$chr[i]),paste(trans.data$chr[i])),
                              res=10000,
                              search.size=search.size)
    
    
      averages[i]=RS$confInterval[2]
    
      totals=as.vector(na.omit(RS$resampled.totals))
      
      vec=ifelse(totals>reads[i], 1, 0)
    
      p.values[i]=sum(na.omit(vec))/length(totals)
      
    }
    
    
    
    
    
  }
  
  #summarize
  info.list = cbind.data.frame(trio.attr$Attributes$cis$trio.idx, 
                               reads, averages, p.values,
                               trio.attr$Attributes$cis$chr,
                               trio.attr$Attributes$trans$chr,
                               trio.attr$Attributes$cis$variant_pos,
                               trio.attr$Attributes$trans$left,
                               trio.attr$Attributes$trans$right,
                               hic.extent)
  
  colnames(info.list)=c("trio.idx", "reads", "expected", "P(>obs)", "cis.chr", 
                        "trans.chr","variant.pos", "trans.left","trans.right", 
                        "variant.lower.bound", "variant.upper.bound",
                        "trans.lower.bound", "trans.upper.bound")
  
  return(info.list)
  
  
  
}

#==================================END_Function========================================================





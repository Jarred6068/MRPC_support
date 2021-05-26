
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
  #A simple helper function used by get_trio_attr() to find trans gene attributes
  
  tissue.names=read.csv("/mnt/ceph/jarredk/AddisReRunFiles/tissuenames.csv", header = T)

    meta2=read.delim(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissue.names[1,2],
                    "_AllPC/",tissue.names[1,2],".genepos.V8.txt",sep="")) 
      
      TGattr=subset(meta2, gene_id==TransGeneName)
  
  return(TGattr)
  
}

#==============================================Helper_Function==============================================

find_cis=function(CisGeneName=NULL, tiss.name=NULL){
  #A simple helper function used by get_trio_attr() to find trans gene attributes
  name=subset(tissue.names, tissue.name1==tiss.name)
  meta=loadRData(fileName = paste("/mnt/ceph/jarredk/AddisReRunFiles/", name[3], ".eGenes.V8.unique.RData", sep = ""))
  
  CGattr=subset(meta, gene_id==CisGeneName)
  
  return(CGattr)
  
}



#==================================Begin_Helper_Function_2==================================================
get_trio_attr=function(trio.index=NULL, tissue.name="CellsEBVtransformedlymphocytes", verbose=TRUE){
  #This function searches for the specified a trios attributes 
  #-------------------------------------------------------------------------------------------------
  #SYNTAX: 
  #trio.index -- the index of the trio for which attributes are desired
  #tissue.name -- the tissue of the trio 
  
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
                               verbose=FALSE, plot.hist=FALSE, trio=NULL, tiss=NULL, FDR=NULL){
  
  #SYNTAX:
  #filepath - feeds to extrac_hic to obtain min/max bounds for both chromosomes
  #chrs - the number, or X/Y identifying the chromosome of variant and trans gene
  #res - the desired resolution of the extraction
  #search.size - specifies the size of bin around chrs position to be searched
  #resamples - the number of resampling iterations
  #verbose - logical indicating if bounds should be printed
  #plot.hist - plots and saves to file the histogram of resampled interactions:
  
  #------if--plot.hist--TRUE: the following must be specified
  #trio - the number identifying the trio of interest: i.e 3324
  #tiss - the tissue for the current trio and interaction data: i.e "SkingNotSunExposed"
  #FDR - the method of FDR in MRPC by which the trio was classified:"LOND" or "ADDIS"
  
  
  #for histogram plotting
  tissues=c("CellsCulturedfibroblasts","SkinNotSunExposed","Lung","CellsEBVtransformedlymphocytes")
  tisspath=c("fibroblast_cells","Skin","Lung","lymphoblastoid_cells")
  
  reads=rep(0, resamples)
  sampler=NULL
  CI=NULL
  i=1
  data.hic2=extract_hic(fileName = filePath, 
                        chrs = chrs,
                        resol=res)
  #obtain the bounds of the chromosomes
  bounds=c(min(data.hic2$x), max(data.hic2$x), min(data.hic2$y), max(data.hic2$y))
  
  if(verbose==TRUE){
    print(bounds)
  }

  for( i in 1:length(reads)){
    #randomly select a position on the chromosomes
    gene1=runif(2, min = bounds[1]+search.size, max = bounds[2]-search.size)
    gene2=runif(2, min = bounds[3]+search.size, max = bounds[4]-search.size)
    
    #establish snp - gene search "box"
    left.pos=paste(chrs[1],":",
                   gene1[1]-search.size,":",
                   gene1[1]+search.size, sep = "")
  
    right.pos=paste(chrs[2],":",
                    gene2[2]-search.size,":",
                    gene2[2]+search.size, sep = "")
  
    #extract reads in search box for each gene/variant
    hic.obj=extract_hic(fileName=filePath, chrs = c(left.pos, right.pos), resol = res)
    
    #count interactions
    reads[i]=ifelse(empty(as.data.frame(hic.obj))==TRUE, NA, sum(hic.obj$counts))
    
  }
  
  #histogram 
  if(plot.hist==TRUE){
    
    idx=which(tiss==tissues)
    G=tisspath[idx]
    
    png(paste("/mnt/ceph/jarredk/HiC_Analyses/",FDR,"/Histograms/",G ,"/rplot_", trio ,".png", sep = ""))
    
    H1=hist(reads, breaks = 10)
    plot(H1)
    
    dev.off()
    
  }
  
  reads=na.omit(reads)
  
  #calculate confidence interval around mean number of interactions
  X.bar=mean(reads)
  X.sd=sd(reads)
  lower.limit=X.bar-abs(qnorm(0.025))*X.sd
  upper.limit=X.bar+abs(qnorm(0.025))*X.sd
  CI=c(lower.limit, X.bar, upper.limit)
  names(CI)=c("lower.limit", "mean", "upper.limit")
  
  #return resampled data (and hist object if plot.hist==TRUE)
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
                           verbose=TRUE, pack.path="/mnt/ceph/jarredk/Rpackages", plot.h=FALSE, FDR=NULL){
  #SYNTAX:
  #hic.filename -- the path to the desired .hic file
  #trio -- the trios index (or vector of indexes) typically obtained from ADDIS.M1.check() or ADDIS.PostProc() 
  #resolution -- the BP resolution desired and to be passed to "extract_hic() and Resample_interactions()"
  #search.size -- the size of the bin around gene positions: also passed to Resample_interactions()
  #tiss -- a string specifying the name of the tissue to which the trio belongs 
  #Verbose -- logical which if true prints checkpoints as trios are checked
  #pack.path -- the path to you Rpackages library for dependencies such as StrawR
  #plot.h -- logical passed to Resample_interactions()
  #FDR -- string specifying "LOND" or "ADDIS" passed to Resample_interactions()
  
  
  
  #obtain trio attributes
  trio.attr=get_trio_attr(trio.index=trios, tissue.name = tiss, verbose = verbose)
  
  cis.data=trio.attr$Attributes$cis
  trans.data=trio.attr$Attributes$trans
  #pre-allocate
  idx = length(trios)
  reads = NULL
  averages = NULL
  p.values = NULL
  p.values2 = NULL
  total.nas = NULL
  resampled_dataset=list()
  
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
      #if no reads for variant<-->gene return NA's
      averages[i]=NA
      p.values[i]=NA
      p.values2[i]=NA
    }else{
      #else return the resampled data and determine P(>=obs)
      RS=Resample_interactions(filePath = hic.filename,
                              chrs=c(paste(cis.data$chr[i]),paste(trans.data$chr[i])),
                              res=10000,
                              search.size=search.size,
                              tiss = tiss,
                              trio = trios[i],
                              plot.hist = plot.h,
                              FDR = FDR)
      
      # obtain resampled data, extract NA's and non-NA data points and calc pvalue
      averages[i]=RS$confInterval[2]
      totals=as.vector(na.omit(RS$resampled.totals))
      num_nas=as.vector(attr(RS$resampled.totals,"na.action"))
      
      #MC-integration
      vec=ifelse(totals>=reads[i], 1, 0)
      #P(>obs) = P(>obs and NA) + P(>obs and !NA) = 0 + P(>obs)*P(!NA)
      p.values[i]=( sum(na.omit(vec))/length(totals) )
      p.values2[i]=( sum(na.omit(vec))/(length(totals)+length(num_nas)) )
      
      #return the resampled data set
      resampled_dataset[[i]]=list(sampled=totals, nas=num_nas)
      total.nas[i]=length(num_nas)
      
    }
    
    
    
    
    
  }
  
  
  
  #preallocate space for FDR/FWER adjustments
  bh.thresh=rep(0, length(p.values))
  qvals=rep(0, length(p.values))
  BY=rep(0, length(p.values))
  
  bh.thresh2=rep(0, length(p.values2))
  qvals2=rep(0, length(p.values2))
  BY2=rep(0, length(p.values2))
  
  print(dim(hic.extent))
  print("made it this far1")
  
  #summarize
  info.list = cbind.data.frame(trio.attr$SNPS,
                               trio.attr$Attributes$cis$gene_name,
                               trio.attr$Attributes$cis$gene_id,
                               reads, averages, total.nas, p.values, BY,
                               bh.thresh, qvals, trio.attr$Attributes$cis$chr,
                               trio.attr$Attributes$trans$gene_id,
                               trio.attr$Attributes$trans$chr,
                               trio.attr$Attributes$cis$variant_pos,
                               trio.attr$Attributes$trans$left,
                               trio.attr$Attributes$trans$right,
                               hic.extent)
  
  print(info.list)
  
  info.list2 = cbind.data.frame(trio.attr$SNPS, 
                                trio.attr$Attributes$cis$gene_name,
                                trio.attr$Attributes$cis$gene_id,
                                reads, averages, total.nas, p.values2, BY2,
                                bh.thresh2, qvals2, trio.attr$Attributes$cis$chr,
                                trio.attr$Attributes$trans$gene_id,
                                trio.attr$Attributes$trans$chr,
                                trio.attr$Attributes$cis$variant_pos,
                                trio.attr$Attributes$trans$left,
                                trio.attr$Attributes$trans$right,
                                hic.extent)
  #name cols
  cnames=c("SNP", "gene_name", "cis.gene.ID","obs.reads", "expected", "total_NA's", "P(>obs)", "BY","HB.Adjusted", "qvals",
                        "cis.chr", "trans.gene.ID", "trans.chr","variant.pos", "trans.left","trans.right", "variant.lower.bound", 
                        "variant.upper.bound", "trans.lower.bound", "trans.upper.bound")
  
  colnames(info.list)=cnames
  colnames(info.list2)=cnames
  
  
  #Holm-Bonferroni correction at FWER alpha = 0.05
  #initialize
  alpha=0.05
  m=length(na.omit(p.values))
  sorted.p=sort(p.values, index.return=TRUE, na.last = TRUE)
  sorted.p2=sort(p.values2, index.return=TRUE, na.last = TRUE)
  print(sorted.p$x)
  HB.adjust=rep(0, length(p.values))
  
  #calculate the rejections using step-down procedure
  for(k in 1:m){
    HB.adjust[k]=sorted.p$x[k]*(m+1-k)
  }
  
  #store in output table
  info.list=info.list[sorted.p$ix,]
  info.list2=info.list2[sorted.p2$ix,]
  HB.adjust=ifelse(HB.adjust>1, 1, HB.adjust)
  HB.adjust=ifelse(HB.adjust==0, NA, HB.adjust)
  info.list$HB.Adjusted=HB.adjust
  
  info.list2$HB.Adjusted=p.adjust(sorted.p2$x, method = "holm")
  
   #get qvalues using BH step-up procedure
  library(qvalue, lib=pack.path)
  
  Q=qvalue(as.vector(na.omit(info.list$`P(>obs)`)), fdr.level = alpha, pi0 = 1)
  Q2=qvalue(as.vector(na.omit(info.list2$`P(>obs)`)), fdr.level = alpha, pi0 = 1)
  
  Qq=as.vector(Q$qvalues)
  Qq2=as.vector(Q2$qvalues)
  
  Qq=ifelse(is.na(info.list$`P(>obs)`)==TRUE, NA, Qq)
  Qq2=ifelse(is.na(info.list2$`P(>obs)`)==TRUE, NA, Qq2)
  info.list$qvals=Qq
  info.list2$qvals=Qq2
  
  #include the Hochberg step up correction
  
  info.list$BY=p.adjust(sorted.p$x, method="BY")
  info.list2$BY=p.adjust(sorted.p2$x, method = "BY")
  
  #return as list
  return(list(summary.table1=info.list, summary.table2=info.list2, data=resampled_dataset))
  
  
  
}

#==================================END_Function========================================================





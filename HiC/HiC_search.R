
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
      
    print("made it this far!")
    #get cis-trans gene attributes
    TGattr=find_trans(TransGeneName=gene.name[2])
    attri.mat.cis[i,] = subset(meta, gene_id == gene.name[1])
    attri.mat.trans[i,] = TGattr
      
  }
  
  
  colnames(attri.mat.trans)=colnames(meta2)
  colnames(attri.mat.cis)=colnames(meta)
  data.list=list(trans=attri.mat.trans, cis=attri.mat.cis)
    
  #check if snp matches
  if(verbose==TRUE){
    ifelse(all.equal(snps, data.list$trans)==TRUE, "ERROR", "Checks.out")
  }
  
  
  
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









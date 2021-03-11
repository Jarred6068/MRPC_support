

#===========================================================================================================

#A function to extract all associated tissues and M classes for a given trio


uncover_net=function(trio_number=NULL, FDR="ADDIS", verbose=TRUE){
  #SYNTAX:
  #trio_number -- the index of the trio you wish to search
  #FDR -- designates which list data to search in (LOND or ADDIS)
  
  
  #load in previous files
  source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")
  
  tissues.vec=tissue.names[,2]
  
  data1=as.data.frame(matrix(0, nrow=length(tissues.vec), ncol = 2))
  names(data1)=c("Tissue", "Model")
  logivec=NULL
  not_in_tiss=rep(FALSE, 5)
  
  
  for( t in 1:length(tissues.vec)){
    
    trioData=loadRData(fileName = paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/", 
                                         tissues.vec[t],"_AllPC/data.snp.cis.trans.final.",tissues.vec[t],".V8.unique.snps.RData", sep = ""))
    
    addisData=ADDIS.M1.check(tissue.names=tissues.vec[t])
  
    #get the M-types (M0, M1, ...)
    Mtypes=names(L1.lond)
    
    for(i in 1:length(L1.lond)){
      
      #match the trio number to trios in each list under Mtypes
      if(FDR=="ADDIS"){

        matched=match(trio_number, L2.addis[[i]])
        
        logivec[i]=ifelse(is.na(matched)==TRUE, FALSE, TRUE)
        
        
      }else{
        
        matched=match(trio_number, L1.lond[[i]])
        
        logivec[i]=ifelse(is.na(matched)==TRUE, FALSE, TRUE)
        
        
      }
      
      
    }
    
    #extract tissues the trio was found in and its Mtype for each tissue
    indicator=isTRUE(all.equal(logivec, not_in_tiss, check.attributes=FALSE))
    
    if(verbose==TRUE & indicator==FALSE){
      print(paste("Found a Match in", tissues.vec[t], "!", sep = " "))
      }
    
    data1[t,1]=ifelse(isTRUE(indicator)==TRUE, NA, tissues.vec[t])
    
    data1[t,2]=ifelse(is.na(match(TRUE, logivec))==TRUE, NA, Mtypes[match(TRUE,logivec)] )
    

  }
  
  #printing
  if(verbose==TRUE){
    print(summary(as.factor(data1$Model)))
    print(paste("Trio found in ",length(na.omit(data1$Tissue)), " Tissues", sep = ""))
    
  }
  
  
  return(data1)
  
}


#===========================================================================================================








































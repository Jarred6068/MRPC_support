

#===========================================================================================================

#A function to extract all associated tissues and M classes for a given trio


uncover_net=function(trio_number=NULL, FDR="ADDIS", verbose=TRUE){
  #SYNTAX:
  #trio_number -- 
  #FDR -- 
  
  
  #load in previous files
  source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")
  
  tissues.vec=tissue.names[,2]
  
  data1=as.data.frame(matrix(0, nrow=length(tissues.vec), ncol = 2))
  names(data1)=c("Tissue", "Model")
  logivec=NULL
  not_in_tiss=rep(FALSE, 5)
  
  
  for( t in 1:length(tissues.vec)){
    
    #Load needed Files
    L1.lond=loadRData(fileName=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/", tissues.vec[t], 
                                       "_AllPC/List.models.", tissues.vec[t], ".all.RData", sep = ""))
    
    L2.addis=loadRData(fileName=paste("/mnt/ceph/jarredk/AddisReRunFiles/List.models.",
                                        tissues.vec[t], ".all.RData", sep = ""))
    
    Mtypes=names(L1.lond)
    
    for(i in 1:length(L1.lond)){
      
      
      if(FDR=="ADDIS"){

        matched=match(trio_number, L2.addis[[i]])
        
        logivec[i]=ifelse(is.na(matched)==TRUE, FALSE, TRUE)
        
        
      }else{
        
        matched=match(trio_number, L1.lond[[i]])
        
        logivec[i]=ifelse(is.na(matched)==TRUE, FALSE, TRUE)
        
        
      }
      
      
    }
    
    indicator=isTRUE(all.equal(logivec, not_in_tiss, check.attributes=FALSE))
    
    if(verbose==TRUE & indicator==FALSE){
      print(paste("Found a Match in", tissues.vec[t], "!", sep = " "))
      }
    
    data1[t,1]=ifelse(isTRUE(indicator)==TRUE, NA, tissues.vec[t])
    
    data1[t,2]=ifelse(is.na(match(TRUE, logivec))==TRUE, NA, Mtypes[match(TRUE,logivec)] )
    

  }
  
  if(verbose==TRUE){
    print(summary(as.factor(data1$Model)))
    print(paste("Trio found in ",length(na.omit(data1$Tissue)), " Tissues", sep = ""))
  }
  return(data1)
  
  
}











































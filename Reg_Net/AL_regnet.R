#load in previous files
# source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")
# 
# tissues.vec=tissue.names[,2]
# for(t in 1:length(tissues.vec)){
#   M1.Data=ADDIS.M1.check(tissue.names=tissues.vec[t])
#   
#   save(M1.Data, file = paste("/mnt/ceph/jarredk/Reg_Net/AL_M1Data/",tissues.vec[t],".Rdata",sep = ""))
# 
# }


#===========================================================================================================

#A function to extract all associated tissues and M classes for a given trio


uncover_net=function(tissue.selection=1, FDR="ADDIS", mediator="trans", verbose=TRUE){
  #SYNTAX:
  #trio_number -- the index of the trio you wish to search
  #FDR -- designates which list data to search in (LOND or ADDIS)
  
  
  #load in previous files
  source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")
  
  tissues.vec=tissue.names[tissue.selection,2]
  all.tissues=tissue.names[,2]
  loaded.datasets=list()
  
  #loading bar
  lbar <- txtProgressBar(min = 0, max = length(all.tissues), style = 3)
  #load in all datasets
  for(i in 1:length(all.tissues)){
    
    loaded.datasets[[i]]=loadRData(fileName = paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/", 
                                                    all.tissues[i],"_AllPC/data.snp.cis.trans.final.",all.tissues[i],".V8.unique.snps.RData", sep = ""))
    
    Sys.sleep(0.05)
    setTxtProgressBar(lbar, i)
    
  }
  #close loading bar
  close(lbar)
  print("--Data:--loading--complete--")
  
  #preallocate space
  matchtab1=as.data.frame(matrix(0, nrow = length(tissues.vec), ncol = 2))
  colnames(matchtab1)=c("tissue", "found")
  all.matches=list()
  
  #begin for all tissues in tissues selection
  for( t in 1:length(tissues.vec)){
    
    trioData=loadRData(fileName = paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/", 
                                         tissues.vec[t],"_AllPC/data.snp.cis.trans.final.",tissues.vec[t],".V8.unique.snps.RData", sep = ""))
    
    M1.Data=loadRData(fileName = paste("/mnt/ceph/jarredk/Reg_Net/AL_M1Data/", tissues.vec[t],".Rdata", sep = ""))
    
    addisM1 = M1.Data$Catalog[[1]][[1]][[2]]
    londM1 =  M1.Data$Catalog[[1]][[1]][[1]]
#=========================================================
#========================ADDIS============================
    if(FDR=="ADDIS"){
      
#-------------------Addis-Trans---------------------------
      if(mediator=="trans"){
        
        matchtab2=list()
        
        for(j in 1:length(addisM1$type2)){
          
          index=(addisM1$type2[j]-1)*3+(1:3)
          members=colnames(trioData[,index])[3]
          print(members)
          
          for(k in 1:length(all.tissues)){
            
            trioData2=colnames(loaded.datasets[[k]])
            
            matched=match(members, trioData2)
            #print(matched)
            
            matchtab1[k,1]=ifelse(is.na(matched)==TRUE, NA, all.tissues[k])
            matchtab1[k,2]=ifelse(is.na(matched)==TRUE, 0, 1)

          }
          
          matchtab2[[j]]=matchtab1
          #print(matchtab2)
        }
        
        all.matches[[t]]=matchtab2
        #print(all.matches)
        
#------------------addis-cis------------------------        
      }else{
        
        matchtab2=list()
        
        for(j in 1:length(addisM1$type1)){
          
          index=(addisM1$type1[j]-1)*3+(1:3)
          members=colnames(trioData[,index])[3]
          print(members)
          
          for(k in 1:length(all.tissues)){
            
            trioData2=colnames(loaded.datasets[[k]])
            
            matched=match(members, trioData2)
            #print(matched)
            
            matchtab1[k,1]=ifelse(is.na(matched)==TRUE, NA, all.tissues[k])
            matchtab1[k,2]=ifelse(is.na(matched)==TRUE, 0, 1)
            
          }
          
          matchtab2[[j]]=matchtab1
          #print(matchtab2)
        }
        
        all.matches[[t]]=matchtab2
        #print(all.matches)
        
      }
      
#=====================================================
#======================LOND===========================

    }else{
      
#------------------lond---trans-----------------------   
      if(mediator=="trans"){
        
        matchtab2=list()
        
        for(j in 1:length(londM1$type2)){
          
          index=(londM1$type2[j]-1)*3+(1:3)
          members=colnames(trioData[,index])[3]
          print(members)
          
          for(k in 1:length(all.tissues)){
            
            trioData2=colnames(loaded.datasets[[k]])
            
            matched=match(members, trioData2)
            #print(matched)
            
            matchtab1[k,1]=ifelse(is.na(matched)==TRUE, NA, all.tissues[k])
            matchtab1[k,2]=ifelse(is.na(matched)==TRUE, 0, 1)
            
          }
          
          matchtab2[[j]]=matchtab1
          #print(matchtab2)
        }
        
        all.matches[[t]]=matchtab2
        #print(all.matches)
      
     
#-------------------lond--cis------------------------       
      }else{
        
        matchtab2=list()
        
        for(j in 1:length(londM1$type1)){
          
          index=(londM1$type1[j]-1)*3+(1:3)
          members=colnames(trioData[,index])[3]
          print(members)
          
          for(k in 1:length(all.tissues)){
            
            trioData2=colnames(loaded.datasets[[k]])
            
            matched=match(members, trioData2)
            #print(matched)
            
            matchtab1[k,1]=ifelse(is.na(matched)==TRUE, NA, all.tissues[k])
            matchtab1[k,2]=ifelse(is.na(matched)==TRUE, 0, 1)
            
          }
          
          matchtab2[[j]]=matchtab1
          #print(matchtab2)
        }
        
        all.matches[[t]]=matchtab2
        #print(all.matches)
        
      }
      
    }

  } #closes first 4loop
  
  names(all.matches)=tissues.vec
  
  if(FDR=="ADDIS" & mediator == "trans"){names(all.matches[[1]])=as.character(addisM1$type2)}
  if(FDR=="ADDIS" & mediator == "cis"){names(all.matches[[1]])=as.character(addisM1$type1)}
  if(FDR=="LOND" & mediator == "trans"){names(all.matches[[1]])=as.character(londM1$type2)}
  if(FDR=="LOND" & mediator == "cis"){names(all.matches[[1]])=as.character(londM1$type1)}
  
  return(all.matches)
  
} #closes func


#===========================================================================================================




#source("/mnt/ceph/jarredk/Reg_Net/AL_regnet.R")



































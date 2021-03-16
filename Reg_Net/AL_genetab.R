#save all M1 data for speed
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


sort.gt=function(tissue.selection=1, FDR="ADDIS", mediator="trans", save.data=FALSE, verbose=TRUE){
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
  all.matches=vector(mode = "list", length = length(tissues.vec))
  all.tables=vector(mode = "list", length = length(tissues.vec))
  #pre-naming
  colnames(matchtab1)=c("tissue", "found")
  names(all.matches)=tissues.vec
  names(all.tables)=tissues.vec
  
  #loading bar 2
  print("---Running:--sort.gt()---")
  lbar2 <- txtProgressBar(min = 0, max = length(tissues.vec), style = 3)
  #begin for all tissues in tissues selection
  for( t in 1:length(tissues.vec)){
    
    Sys.sleep(0.005)
    setTxtProgressBar(lbar2, t)
    #read in the triodata for selected
    trioData=loadRData(fileName = paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/", 
                                         tissues.vec[t],"_AllPC/data.snp.cis.trans.final.",tissues.vec[t],".V8.unique.snps.RData", sep = ""))
    
    #read in M1-types data from ADDIS.M1.check() 
    M1.Data=loadRData(fileName = paste("/mnt/ceph/jarredk/Reg_Net/AL_M1Data/", tissues.vec[t],".Rdata", sep = ""))
    #get M1's
    addisM1 = M1.Data$Catalog[[1]][[1]][[2]]
    londM1 =  M1.Data$Catalog[[1]][[1]][[1]]
#=========================================================
#========================ADDIS============================
    if(FDR=="ADDIS"){
      
#-------------------Addis-Trans---------------------------
      if(mediator=="trans"){
        
        matchtab2=vector(mode = "list", length = length(all.tissues))
        names(matchtab2)=all.tissues
        summary.table1=as.data.frame(matrix(0, nrow = length(addisM1$type2), ncol = 3))
        colnames(summary.table1)=c("Count.Tissues", "Percent.Tissues","Gene.Name")
        gene.name=NULL
        
        for(j in 1:length(addisM1$type2)){
          
          index=(addisM1$type2[j]-1)*3+(1:3)
          members=colnames(trioData[,index])[3]
          if(verbose==TRUE){print(members)}
          gene.name[j]=members
          
          for(k in 1:length(all.tissues)){
            
            trioData2=colnames(loaded.datasets[[k]])
            
            matched=match(members, trioData2)
            #print(matched)
            
            matchtab1[k,1]=ifelse(is.na(matched)==TRUE, NA, all.tissues[k])
            matchtab1[k,2]=ifelse(is.na(matched)==TRUE, 0, 1)
            matchtab2[[k]][j]=ifelse(is.na(matched)==TRUE, NA, members)

          }
          
          summary.table1[j, 1]=sum(matchtab1[,2])
          summary.table1[j, 2]=sum(matchtab1[,2])/length(all.tissues)
          summary.table1[j, 3]=gene.name[j]

        }

        all.matches[[t]]=matchtab2
        all.tables[[t]]=summary.table1
        
#------------------addis-cis------------------------        
      }else{
        
        matchtab2=vector(mode = "list", length = length(all.tissues))
        names(matchtab2)=all.tissues
        summary.table1=as.data.frame(matrix(0, nrow = length(addisM1$type1), ncol = 3))
        colnames(summary.table1)=c("Count.Tissues", "Percent.Tissues","Gene.Name")
        gene.name=NULL
        
        for(j in 1:length(addisM1$type1)){
          
          index=(addisM1$type1[j]-1)*3+(1:3)
          members=colnames(trioData[,index])[3]
          if(verbose==TRUE){print(members)}
          gene.name[j]=members
          
          for(k in 1:length(all.tissues)){
            
            trioData2=colnames(loaded.datasets[[k]])
            
            matched=match(members, trioData2)
            #print(matched)
            
            matchtab1[k,1]=ifelse(is.na(matched)==TRUE, NA, all.tissues[k])
            matchtab1[k,2]=ifelse(is.na(matched)==TRUE, 0, 1)
            matchtab2[[k]][j]=ifelse(is.na(matched)==TRUE, NA, members)
            
          }
          
          summary.table1[j, 1]=sum(matchtab1[,2])
          summary.table1[j, 2]=sum(matchtab1[,2])/length(all.tissues)
          summary.table1[j, 3]=gene.name[j]
          
        }
        
        all.matches[[t]]=matchtab2
        all.tables[[t]]=summary.table1
        
      }
      
#=====================================================
#======================LOND===========================

    }else{
      
#------------------lond---trans-----------------------   
      if(mediator=="trans"){
        
        matchtab2=vector(mode = "list", length = length(all.tissues))
        names(matchtab2)=all.tissues
        summary.table1=as.data.frame(matrix(0, nrow = length(londM1$type2), ncol = 3))
        colnames(summary.table1)=c("Count.Tissues", "Percent.Tissues","Gene.Name")
        gene.name=NULL
        
        for(j in 1:length(londM1$type2)){
          
          index=(londM1$type2[j]-1)*3+(1:3)
          members=colnames(trioData[,index])[3]
          if(verbose==TRUE){print(members)}
          gene.name[j]=members
          
          for(k in 1:length(all.tissues)){
            
            trioData2=colnames(loaded.datasets[[k]])
            
            matched=match(members, trioData2)
            #print(matched)
            
            matchtab1[k,1]=ifelse(is.na(matched)==TRUE, NA, all.tissues[k])
            matchtab1[k,2]=ifelse(is.na(matched)==TRUE, 0, 1)
            matchtab2[[k]][j]=ifelse(is.na(matched)==TRUE, NA, members)
            
          }
          
          summary.table1[j, 1]=sum(matchtab1[,2])
          summary.table1[j, 2]=sum(matchtab1[,2])/length(all.tissues)
          summary.table1[j, 3]=gene.name[j]
          
        }

        all.matches[[t]]=matchtab2
        all.tables[[t]]=summary.table1
      
     
#-------------------lond--cis------------------------       
      }else{
        
        matchtab2=vector(mode = "list", length = length(all.tissues))
        names(matchtab2)=all.tissues
        summary.table1=as.data.frame(matrix(0, nrow = length(londM1$type1), ncol = 3))
        colnames(summary.table1)=c("Count.Tissues", "Percent.Tissues","Gene.Name")
        gene.name=NULL
        
        for(j in 1:length(londM1$type1)){
          
          index=(londM1$type1[j]-1)*3+(1:3)
          members=colnames(trioData[,index])[3]
          if(verbose==TRUE){print(members)}
          gene.name[j]=members
          
          for(k in 1:length(all.tissues)){
            
            trioData2=colnames(loaded.datasets[[k]])
            
            matched=match(members, trioData2)
            #print(matched)
            
            matchtab1[k,1]=ifelse(is.na(matched)==TRUE, NA, all.tissues[k])
            matchtab1[k,2]=ifelse(is.na(matched)==TRUE, 0, 1)
            matchtab2[[k]][j]=ifelse(is.na(matched)==TRUE, NA, members)
            
          }
          
          summary.table1[j, 1]=sum(matchtab1[,2])
          summary.table1[j, 2]=sum(matchtab1[,2])/length(all.tissues)
          summary.table1[j, 3]=gene.name[j]
          
        }
        
        all.matches[[t]]=matchtab2
        all.tables[[t]]=summary.table1
        
      }
      
    }

  } #closes first 4loop
  
  #close loading bar
  close(lbar2)
  print("--Running:--sort.gt()--complete--")
  
  if(save.data==TRUE){
    print("--Saving-Files-To:--/mnt/ceph/jarredk/Reg_Net/--")
    save.set=list(all.matches, all.tables)
    save(save.set, file = paste("/mnt/ceph/jarredk/Reg_Net/",FDR,"_", mediator,"_", "sorted_data_list.Rdata", sep = ""))
  }
  
  return(list( match.tables=all.matches, summary.tables=all.tables))
  
} #closes func


#===========================================================================================================




#source("/mnt/ceph/jarredk/Reg_Net/AL_genetab.R")
#sort.gt(tissue.selection=c(1:48), FDR="ADDIS", mediator="trans", save.data=TRUE, verbose=FALSE)


#===========================================================================================================

table.create=function(tissues=c(1:48), FDR.method="ADDIS", mediator.type="trans", run.sort.gt=FALSE, verbose=FALSE){
  
  source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")
  
  
  if(run.sort.gt==TRUE){
    
    dat=sort.gt(tissue.selection = tissues, FDR = FDR.method, mediator = mediator.type, verbose = verbose)
    
  }else{
    
    dat=loadRData(fileName=paste("/mnt/ceph/jarredk/Reg_Net/",FDR.method,"_", mediator.type,"_", "sorted_data_list.Rdata", sep = ""))
    
  }
  
  
  prelist=vector(mode = "list", length = length(tissues))
  names(prelist)=tissue.names[tissues,2]
  shared=as.data.frame(matrix(0, nrow = length(tissues), ncol=2))
  colnames(shared)=c("# genes", "tissue")
  nam.vec=tissue.names[tissues,2]
  
  
  #loading bar
  lbar = txtProgressBar(min = 0, max = length(tissues), style = 3)
  
  for(i in 1:length(tissues)){
    
    for(j in 1:length(tissues)){
      
      prelist[[j]]=unique(c(prelist[[j]], as.vector(na.omit(dat[[1]][[i]][[j]]))))
      shared[j,1]=length(prelist[[j]])
      shared[j,2]=nam.vec[j]
      
    }
    
    Sys.sleep(0.05)
    setTxtProgressBar(lbar, i)
    
  }
  
  #close loading bar
  close(lbar)
  print("--Data:--Processing--Stage_1--complete--")
  
  
  return(list(all.genes=prelist, shared.tissues=shared))
  
}































#source("/mnt/ceph/jarredk/Reg_Net/AL_genetabV2.R")

#=============================================================================================================

sort.gt2=function(tissue.selection=1, FDR="ADDIS", mediator="trans", save.data=FALSE, verbose=FALSE){
  #SYNTAX:
  #trio_number -- the index of the trio you wish to search
  #FDR -- designates which list data to search in (LOND or ADDIS)
  
  
  #load in previous files
  source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")
  
  tissues.vec=tissue.names[tissue.selection,2]
  all.tissues=tissue.names[,2]
  loaded.datasets=list()
  
  # #loading bar
  # lbar <- txtProgressBar(min = 0, max = length(all.tissues), style = 3)
  # #load in all datasets
  # for(i in 1:length(all.tissues)){
  #   
  #   loaded.datasets[[i]]=loadRData(fileName = paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/", 
  #                                                   all.tissues[i],"_AllPC/data.snp.cis.trans.final.",all.tissues[i],".V8.unique.snps.RData", sep = ""))
  #   
  #   Sys.sleep(0.05)
  #   setTxtProgressBar(lbar, i)
  #   
  # }
  # #close loading bar
  # close(lbar)
  # print("--Data:--loading--complete--")
  # 
  #load in meta data from BIOMart
  meta.data=read.table(file="/mnt/ceph/jarredk/Reg_Net/mart_export.txt", sep=",", header=T)
  
  #loading bar 2
  print("---Running:--sort.gt2()---")
  lbar2 <- txtProgressBar(min = 0, max = length(tissues.vec), style = 3)
  
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
    
    #allocate space:
    master=vector(mode = "list", length = length(all.tissues))
    names(master)=all.tissues
    
    #=========================================================
    #========================ADDIS============================
    if(FDR=="ADDIS"){
      
      #-------------------Addis-Trans---------------------------
      if(mediator=="trans"){
        trio.names=as.data.frame(matrix(0, nrow = length(addisM1$type2), ncol = 9))
        colnames(trio.names)=c("SNP", "Cis.Gene.ID","Trans.Gene.ID", "Mediator", colnames(meta.data)[2:6])
        
        for(i in 1:length(addisM1$type2)){
          
          index=(addisM1$type2[i]-1)*3+(1:3)
          mediator.name=colnames(trioData[,index])[3]
          
          mediator.meta=meta.data[which(meta.data$Gene.stable.ID==mediator.name),2:6]
          
          trio.names[i,]=c(colnames(trioData[,index]),"M1T2", mediator.meta)
          
          
          
        }
        
        master[[t]]=trio.names
        print(trio.names)

        
        #------------------addis-cis------------------------        
      }else{
        
        trio.names=as.data.frame(matrix(0, nrow = length(addisM1$type1), ncol = 9))
        colnames(trio.names)=c("SNP", "Cis.Gene.ID","Trans.Gene.ID", "Mediator", colnames(meta.data)[2:6])
        
        for(i in 1:length(addisM1$type1)){
          
          index=(addisM1$type2[i]-1)*3+(1:3)
          mediator.name=colnames(trioData[,index])[3]
          
          mediator.meta=meta.data[which(meta.data$Gene.stable.ID==mediator.name),2:6]
          
          trio.names[i,]=c(colnames(trioData[,index]),"M1T2", mediator.meta)
          
          
          
        }
        
        master[[t]]=trio.names
        print(trio.names)
        
        
        
      }
      
      #=====================================================
      #======================LOND===========================
      
    }else{
      
      #------------------lond---trans-----------------------   
      if(mediator=="trans"){
        
          

        
        
        #-------------------lond--cis------------------------       
      }else{
        
        
        
      }
      
    }
    
  } #closes first 4loop
  
  #close loading bar
  close(lbar2)
  print("--Running:--sort.gt()--complete--")
  
  
  return(master)
  
  
  
}
  
  
  #=============================================================================================================
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
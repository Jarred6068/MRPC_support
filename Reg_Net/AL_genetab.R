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

#===================================================================================================
ADDIS.M1.checkV2=function(trio.index=NULL, tissue.name=NULL, FDR="ADDIS", verbose=FALSE){
  #SYNTAX: 
  #tissue.names -- can be either a single tissue name or character vector of names
  #corresponding to the format of column 2 in tissuenames.csv
  
  typeM=NULL
  

  
  #establish ground truths for M1.1 and M1.2
  Truth.M1.type1=matrix(c(1,0,0,1,0,0), nrow = 3, ncol = 2, byrow = T)
  Truth.M1.type2=matrix(c(0,1,0,0,1,0), nrow = 3, ncol = 2, byrow = T)
  

      if(FDR=="ADDIS"){
        
        results=Lond2Addis.lookup(trio.index = trio.index, tissue.name = tissue.name)
        MRPC.graph=results$MRPC.table.ADDIS[1:3,2:3]
        
        typeM=ifelse(all.equal(MRPC.graph, Truth.M1.type1, check.attributes=FALSE)==TRUE, "M1.1", NA)
        typeM=ifelse(all.equal(MRPC.graph, Truth.M1.type2, check.attributes=FALSE)==TRUE, "M1.2", typeM)
        
      }else{
        
        results=Lond2Addis.lookup(trio.index = trio.index, tissue.name = tissue.name)
        MRPC.graph=results$MRPC.table.LOND[1:3,2:3]
        
        typeM=ifelse(all.equal(MRPC.graph, Truth.M1.type1, check.attributes=FALSE)==TRUE, "M1.1", NA)
        typeM=ifelse(all.equal(MRPC.graph, Truth.M1.type2, check.attributes=FALSE)==TRUE, "M1.2", typeM)
        
      }
  
  return(typeM)
  
}

#==========================================END--FUNCTION_3===================================


#===========================================================================================================

#A function to extract all associated tissues and M classes for a given trio

#source("/mnt/ceph/jarredk/Reg_Net/AL_genetab.R")
#sort.gt(tissue.selection=c(1:48), FDR="ADDIS", mediator="trans", save.data=TRUE, verbose=FALSE)



sort.gt=function(tissue.selection=1, FDR="ADDIS", mediator="trans", save.data=FALSE, verbose=FALSE){
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
  matchtab1=as.data.frame(matrix(0, nrow = length(tissues.vec), ncol = 3))
  all.matches=vector(mode = "list", length = length(tissues.vec))
  all.matches2=vector(mode = "list", length = length(tissues.vec))
  all.tables=vector(mode = "list", length = length(tissues.vec))
  #pre-naming
  colnames(matchtab1)=c("tissue", "found&M1", "Mtype")
  names(all.matches)=tissues.vec
  names(all.matches2)=tissues.vec
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
        matchtab3=list()
        summary.table1=as.data.frame(matrix(0, nrow = length(addisM1$type2), ncol = 3))
        colnames(summary.table1)=c("Count.Tissues", "Percent.Tissues","Gene.Name")
        gene.name=NULL
        #print(length(addisM1$type2))
        
        for(j in 1:length(addisM1$type2)){
          
          # begin.time=Sys.time()
          index=(addisM1$type2[j]-1)*3+(1:3)
          members=colnames(trioData[,index])[3]
          if(verbose==TRUE){print(members)}
          gene.name[j]=members

          for(k in 1:length(all.tissues)){
            
            any.matches=any(which(colnames(loaded.datasets[[k]])==members))
            if(isTRUE(any.matches)==FALSE){
              matched=NA
            }else{
              matched=which(colnames(loaded.datasets[[k]])==members)
            }
            
            
            if(isTRUE(any.matches)==FALSE){
              
              trio.idx=NA
              typeM=NA
              
            }else{
              
              trio.idx=NULL
              for(i in 1:length(matched)){
                trio.idx[i]=ifelse(isTRUE(matched[i]%%3==0), matched[i]/3, (matched[i]+1)/3)
              }

              M1.Data2=loadRData(fileName = paste("/mnt/ceph/jarredk/Reg_Net/AL_M1Data/", all.tissues[k],".Rdata", sep = ""))
              addis2M1 = M1.Data2$Catalog[[1]][[1]][[2]]
              
              matched2=NULL
              for(i in 1:length(trio.idx)){
                matched2[i]=any(which(addis2M1$type2==trio.idx[i]))
              }
              
              typeM=ifelse(isTRUE(any(matched2))==FALSE, NA, "M1.2")
              print(typeM)
              
            }
            
             print(any(matched2))
             print(matched2)
             print("=================================================")

            matchtab1[k,1]=ifelse(is.na(typeM)==TRUE, NA, all.tissues[k])
            matchtab1[k,2]=ifelse(is.na(typeM)==TRUE, 0, 1)
            matchtab1[k,3]=typeM
            matchtab2[[k]][j]=ifelse(is.na(typeM)==TRUE, NA, members)
            

          }
          
          #print(matchtab1)
          matchtab3[[j]]=matchtab1[,2]
          summary.table1[j, 1]=sum(matchtab1[,2])
          summary.table1[j, 2]=sum(matchtab1[,2])/length(all.tissues)
          summary.table1[j, 3]=gene.name[j]
          
          # end.time=Sys.time()
          # print(end.time-begin.time)

        }

        all.matches[[t]]=matchtab2
        all.matches2[[t]]=matchtab3
        all.tables[[t]]=summary.table1
        
#------------------addis-cis------------------------        
      }else{
        

        

        
      }
      
#=====================================================
#======================LOND===========================

    }else{
      
#------------------lond---trans-----------------------   
      if(mediator=="trans"){
        
        matchtab2=vector(mode = "list", length = length(all.tissues))
        names(matchtab2)=all.tissues
        matchtab3=list()
        summary.table1=as.data.frame(matrix(0, nrow = length(londM1$type2), ncol = 3))
        colnames(summary.table1)=c("Count.Tissues", "Percent.Tissues","Gene.Name")
        gene.name=NULL
        #print(length(addisM1$type2))
        
        for(j in 1:length(londM1$type2)){
          
          # begin.time=Sys.time()
          index=(londM1$type2[j]-1)*3+(1:3)
          members=colnames(trioData[,index])[3]
          if(verbose==TRUE){print(members)}
          gene.name[j]=members
          
          for(k in 1:length(all.tissues)){
            
            any.matches=any(which(colnames(loaded.datasets[[k]])==members))
            if(isTRUE(any.matches)==FALSE){
              matched=NA
            }else{
              matched=which(colnames(loaded.datasets[[k]])==members)
            }
            
            
            if(isTRUE(any.matches)==FALSE){
              
              trio.idx=NA
              typeM=NA
              
            }else{
              
              trio.idx=NULL
              for(i in 1:length(matched)){
                trio.idx[i]=ifelse(isTRUE(matched[i]%%3==0), matched[i]/3, (matched[i]+1)/3)
              }
              
              M1.Data2=loadRData(fileName = paste("/mnt/ceph/jarredk/Reg_Net/AL_M1Data/", all.tissues[k],".Rdata", sep = ""))
              lond2M1 = M1.Data2$Catalog[[1]][[1]][[1]]
              
              matched2=NULL
              for(i in 1:length(trio.idx)){
                matched2[i]=any(which(lond2M1$type2==trio.idx[i]))
              }
              
              typeM=ifelse(isTRUE(any(matched2))==FALSE, NA, "M1.2")
              
            }
            
            # print(any(matched2))
            # print(matched2)
            # print(typeM)
            # print("=================================================")
            
            matchtab1[k,1]=ifelse(is.na(typeM)==TRUE, NA, all.tissues[k])
            matchtab1[k,2]=ifelse(is.na(typeM)==TRUE, 0, 1)
            matchtab1[k,3]=typeM
            matchtab2[[k]][j]=ifelse(is.na(typeM)==TRUE, NA, members)
            
            
          }
          
          #print(matchtab1)
          matchtab3[[j]]=matchtab1[,2]
          summary.table1[j, 1]=sum(matchtab1[,2])
          summary.table1[j, 2]=sum(matchtab1[,2])/length(all.tissues)
          summary.table1[j, 3]=gene.name[j]
          
          # end.time=Sys.time()
          # print(end.time-begin.time)
          
        }
        
        all.matches[[t]]=matchtab2
        all.matches2[[t]]=matchtab3
        all.tables[[t]]=summary.table1

        

      
     
#-------------------lond--cis------------------------       
      }else{
        

        
      }
      
    }

  } #closes first 4loop
  
  #close loading bar
  close(lbar2)
  print("--Running:--sort.gt()--complete--")
  
  if(save.data==TRUE){
    print("--Saving-Files-To:--/mnt/ceph/jarredk/Reg_Net/--")
    save.set=list(all.matches, all.matches2, all.tables)
    save(save.set, file = paste("/mnt/ceph/jarredk/Reg_Net/",FDR,"_", mediator,"_", "sorted_data_list.Rdata", sep = ""))
  }
  
  return(list( match.tables=all.matches, match.df=all.matches2, summary.tables=all.tables))
  
} #closes func


#===========================================================================================================



#===========================================================================================================

table.create=function(tissues=c(1:48), FDR.method="ADDIS", mediator.type="trans", run.sort.gt=FALSE, verbose=FALSE){
  
  source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")
  
  
  if(run.sort.gt==TRUE){
    
    dat=sort.gt(tissue.selection = tissues, FDR = FDR.method, mediator = mediator.type, verbose = verbose)
    
  }else{
    
    dat=loadRData(fileName=paste("/mnt/ceph/jarredk/Reg_Net/",FDR.method,"_", mediator.type,"_", "sorted_data_list.Rdata", sep = ""))
    
  }
  
  #originals
  prelist=vector(mode = "list", length = length(tissues))
  names(prelist)=tissue.names[tissues,2]
  #reduced by uniques
  postlist=vector(mode = "list", length = length(tissues))
  names(prelist)=tissue.names[tissues,2]
  #number of genes in each tissue
  shared=as.data.frame(matrix(0, nrow = length(tissues), ncol=2))
  colnames(shared)=c("# genes", "tissue")
  
  nam.vec=tissue.names[,2]
  
  
  #loading bar
  lbar = txtProgressBar(min = 0, max = length(tissues), style = 3)
  
  for(i in 1:length(tissues)){
    
    for(j in 1:length(tissues)){
      
      prelist[[j]]=c(prelist[[j]], as.vector(na.omit(dat[[1]][[i]][[j]])))
      postlist[[j]]=unique(c(postlist[[j]], as.vector(na.omit(dat[[1]][[i]][[j]]))))
      shared[j,1]=length(prelist[[j]])
      shared[j,2]=nam.vec[j]
      
    }
    
    Sys.sleep(0.05)
    setTxtProgressBar(lbar, i)
    
  }
  
  #close loading bar
  close(lbar)
  print("--Data:--Processing--Stage_1--complete--")
  
  df.new=c(1:49)
  
  #loading bar
  lbar2 = txtProgressBar(min = 0, max = length(tissues), style = 3)
  
  for(t in 1:length(tissues)){
    
    df=as.data.frame.list(dat[[2]][[t]])
    df2=apply(t(df),2,as.numeric)
    df2=cbind.data.frame(dat[[3]][[t]][,3], df2)
    df.new=rbind(df.new, df2)
    
    Sys.sleep(0.05)
    setTxtProgressBar(lbar2, t)
    
  }
  
  df.new=df.new[2:dim(df.new)[1],]
  gene.vec=unique(unlist(prelist))
  #print(length(gene.vec))
  
  idx=match(gene.vec, df.new[,1])
  
  df.new=df.new[idx,]
  
  row.names(df.new)=as.character(seq(dim(df.new)[1]))
  colnames(df.new)=c("gene.name", nam.vec)
  
  
  dist1=as.data.frame(matrix(0, nrow=length(tissues), ncol = 2))
  colnames(dist1)=c("# tissues", "count")
  
  X=rowSums(df.new[,-1])
  
  for(i in 1:length(tissues)){
    
    dist1[i,1]=i
    dist1[i,2]=length(X[X==i])
    
  }
  
  
  close(lbar2)
  print("--Data:--Processing--Stage_2--complete--")
  
  return(list(all.genes=postlist, shared.tissues=shared, master.table=df.new, dist.table=dist1))
  
}





#===========================================================================================================

table2.create=function(tissues=c(1:48), FDR.method="ADDIS", mediator.type="trans", run.sort.gt=FALSE, verbose=FALSE){
  
  meta.data=read.table(file="/mnt/ceph/jarredk/Reg_Net/mart_export.txt", sep=",", header=T)
  
  TT=table.create(tissues = tissues, FDR.method = FDR.method, mediator.type=mediator.type, run.sort.gt = FALSE)
  
  
  
  
  
  
}






















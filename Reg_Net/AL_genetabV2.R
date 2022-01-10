
#source("/mnt/ceph/jarredk/Reg_Net/AL_genetabV2.R")

#=============================================================================================================

sort.gt2=function(tissue.selection=c(1:48), FDR="ADDIS", mediator="trans", save.data=FALSE, verbose=FALSE){
  #SYNTAX:
  #trio_number -- the index of the trio you wish to search
  #FDR -- designates which list data to search in (LOND or ADDIS)
  
  
  #load in previous files
  source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")
  
  tissues.vec=tissue.names[tissue.selection,2]
  all.tissues=tissue.names[,2]
  
  #allocate space:
  master=vector(mode = "list", length = length(all.tissues))
  names(master)=all.tissues
  
  #load in gene meta data from BIO-Mart
  meta.data=read.table(file="/mnt/ceph/jarredk/Reg_Net/mart_export_merged_lncRNA_fixed.txt", sep="\t", header=T)
  meta.names=c("Gene.start.bp","Gene.end.bp","Gene.type","Gene.name","chr")
  sum.stat.names=c("Mean.cis","SD.cis","Mean.trans","SD.trans", "Cor.SNP.cis", "Cor.SNP.trans", "Cor.cis.trans")
  
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
    londM1 = M1.Data$Catalog[[1]][[1]][[1]]
    
    #=========================================================
    #========================ADDIS============================
    if(FDR=="ADDIS"){
      
      #-------------------Addis-Trans---------------------------
      if(mediator=="trans"){
        trio.names=as.data.frame(matrix(0, nrow = length(addisM1$type2), ncol = 21))
        colnames(trio.names)=c("SNP", "Cis.Gene.ID","Trans.Gene.ID", "Mediator",
                               paste0("trans.",meta.names), paste0("cis.",meta.names), sum.stat.names)
        
        for(i in 1:length(addisM1$type2)){
          
          index=(addisM1$type2[i]-1)*3+(1:3)
          trio=trioData[,index]
          gene.names=colnames(trio)
          mediator.meta=meta.data[which(meta.data$Gene.stable.ID==gene.names[3])[1],1:5]
          nonmed.meta=meta.data[which(meta.data$Gene.stable.ID==gene.names[2])[1],1:5]
          r=cor(trio, use = "pairwise.complete.obs")
          sum.stats=c(mean(trio[,2]), sd(trio[,2]), mean(trio[,3]), sd(trio[,3]), r[1,2:3], r[2,3])
          trio.names[i,]=c(colnames(trioData[,index]),"AM1T2", mediator.meta, nonmed.meta, sum.stats)
 
        }
        
        master[[t]]=trio.names

        #------------------addis-cis------------------------        
      }else{
        
        trio.names=as.data.frame(matrix(0, nrow = length(addisM1$type1), ncol = 21))
        colnames(trio.names)=c("SNP", "Cis.Gene.ID","Trans.Gene.ID", "Mediator",
                               paste0("cis.",meta.names), paste0("trans.",meta.names), sum.stat.names)
        
        for(i in 1:length(addisM1$type1)){
          
          index=(addisM1$type1[i]-1)*3+(1:3)
          trio=trioData[,index]
          gene.names=colnames(trio)
          mediator.meta=meta.data[which(meta.data$Gene.stable.ID==gene.names[2])[1],1:5]
          nonmed.meta=meta.data[which(meta.data$Gene.stable.ID==gene.names[3])[1],1:5]
          r=cor(trio, use = "pairwise.complete.obs")
          sum.stats=c(mean(trio[,2]), sd(trio[,2]), mean(trio[,3]), sd(trio[,3]), r[1,2:3], r[2,3])
          trio.names[i,]=c(colnames(trioData[,index]),"AM1T1", mediator.meta, nonmed.meta, sum.stats)
          
        }
        
        master[[t]]=trio.names
        
        
      }
      
      #=====================================================
      #======================LOND===========================
      
    }else{
      
      #------------------lond---trans-----------------------   
      if(mediator=="trans"){
        
        trio.names=as.data.frame(matrix(0, nrow = length(londM1$type2), ncol = 21))
        colnames(trio.names)=c("SNP", "Cis.Gene.ID","Trans.Gene.ID", "Mediator",
                               paste0("trans.",meta.names), paste0("cis.",meta.names), sum.stat.names)
        
        for(i in 1:length(londM1$type2)){
          
          index=(londM1$type2[i]-1)*3+(1:3)
          trio=trioData[,index]
          gene.names=colnames(trio)
          mediator.meta=meta.data[which(meta.data$Gene.stable.ID==gene.names[3])[1],1:5]
          nonmed.meta=meta.data[which(meta.data$Gene.stable.ID==gene.names[2])[1],1:5]
          r=cor(trio, use = "pairwise.complete.obs")
          sum.stats=c(mean(trio[,2]), sd(trio[,2]), mean(trio[,3]), sd(trio[,3]), r[1,2:3], r[2,3])
          trio.names[i,]=c(colnames(trioData[,index]),"LM1T2", mediator.meta, nonmed.meta, sum.stats)
          
        }
        
        master[[t]]=trio.names
        
        #-------------------lond--cis------------------------       
      }else{
        
        trio.names=as.data.frame(matrix(0, nrow = length(londM1$type1), ncol = 21))
        colnames(trio.names)=c("SNP", "Cis.Gene.ID","Trans.Gene.ID", "Mediator",
                               paste0("cis.",meta.names), paste0("trans.",meta.names), sum.stat.names)
        
        for(i in 1:length(londM1$type1)){
          
          index=(londM1$type1[i]-1)*3+(1:3)
          trio=trioData[,index]
          gene.names=colnames(trio)
          mediator.meta=meta.data[which(meta.data$Gene.stable.ID==gene.names[2])[1],1:5]
          nonmed.meta=meta.data[which(meta.data$Gene.stable.ID==gene.names[3])[1],1:5]
          r=cor(trio, use = "pairwise.complete.obs")
          sum.stats=c(mean(trio[,2]), sd(trio[,2]), mean(trio[,3]), sd(trio[,3]), r[1,2:3], r[2,3])
          trio.names[i,]=c(colnames(trioData[,index]),"LM1T1", mediator.meta, nonmed.meta, sum.stats)
          
        }
        
        master[[t]]=trio.names
        
      }
      
    }
    
  } #closes first 4loop
  
  #close loading bar
  close(lbar2)
  print("--Running:--sort.gt()--complete--")
  
  
  #save each list as .csv
  if(save.data==TRUE){
    for(i in 1:length(tissues.vec)){
      if(FDR=="ADDIS"){
    
        if(mediator=="trans"){
      
          write.csv(master[[i]], file = paste("/mnt/ceph/jarredk/Reg_Net/AM1T2/",tissues.vec[i],".csv", sep = ""))
      
        }else{
        
          write.csv(master[[i]], file = paste("/mnt/ceph/jarredk/Reg_Net/AM1T1/",tissues.vec[i],".csv", sep = ""))
        
        }
    
      }else{
    
        if(mediator=="trans"){
        
          write.csv(master[[i]], file = paste("/mnt/ceph/jarredk/Reg_Net/LM1T2/",tissues.vec[i],".csv", sep = ""))
      
        }else{
        
          write.csv(master[[i]], file = paste("/mnt/ceph/jarredk/Reg_Net/LM1T1/",tissues.vec[i],".csv", sep = ""))
        
        }
      }
    }
  }
  


  return(master)
  
}

#=============================================================================================================
  
count=function(df=NULL, target=NULL){
  
  unique.genes=unique(df[,target])
  counts.vec=rep(0, length(unique.genes))
  
  for(i in 1:length(unique.genes)){
    
    matched=df[which(df[,target]==unique.genes[i]),]
    counts.vec[i]=length(unique(matched$tissue))
    
  }
  
  counts.vec=cbind.data.frame(unique.genes,counts.vec)
  colnames(counts.vec)=c("Gene.stable.ID", "Num.of.Tissues")
  
  return(counts.vec)
  
}

#=============================================================================================================

binit=function(df=NULL, target=NULL){
  
  cts=count(df=df, target = target)
  unq.cts=c(1:48)
  
  binner=rep(0, length(unq.cts))
  
  for(i in 1:length(binner)){
    
    idx=which(cts$Num.of.Tissues==unq.cts[i])
    binner[i]=length(cts$Gene.stable.ID[idx])
    
  }
  
  binner=cbind.data.frame(binner, unq.cts)
  colnames(binner)=c("Number of Genes", "Number of Tissues Shared")
  
  return(list(counts=cts, bins=binner))
  
}


#=============================================================================================================

match.rowN=function(df=NULL, target.vec=NULL){
  
  idx=rep(0, length(target.vec))
  
  for(i in 1:length(target.vec)){
    
    idx[i]=match(target.vec[i], df)
    
  }
  
  return(idx)
  
}


  
#=============================================================================================================

master.create=function(tissue.selection=c(1:48)){
  #SYNTAX:
  #trio_number -- the index of the trio you wish to search
  #FDR -- designates which list data to search in (LOND or ADDIS)
  
  
  #load in previous files
  source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")
  
  
  tissues.vec=tissue.names[tissue.selection,2]
  all.tissues=tissue.names[,2]
  
  df=rep(0, 15)
  
  
  #loading bar
  lbar <- txtProgressBar(min = 0, max = length(all.tissues), style = 3)
  #load in all datasets
  for(i in 1:length(all.tissues)){
    
    dat1=read.csv(file = paste("/mnt/ceph/jarredk/Reg_Net/AM1T1/", all.tissues[i],".csv",sep = ""))
    tiss1=rep(all.tissues[i], dim(dat1)[1])
    dat1$tissue=tiss1
    
    dat2=read.csv(file = paste("/mnt/ceph/jarredk/Reg_Net/AM1T2/", all.tissues[i],".csv",sep = ""))
    tiss2=rep(all.tissues[i], dim(dat2)[1])
    dat2$tissue=tiss2
    
    dat3=read.csv(file = paste("/mnt/ceph/jarredk/Reg_Net/LM1T1/", all.tissues[i],".csv",sep = ""))
    tiss3=rep(all.tissues[i], dim(dat3)[1])
    dat3$tissue=tiss3
    
    dat4=read.csv(file = paste("/mnt/ceph/jarredk/Reg_Net/LM1T2/", all.tissues[i],".csv",sep = ""))
    tiss4=rep(all.tissues[i], dim(dat4)[1])
    dat4$tissue=tiss4
    
    df=rbind.data.frame(df, dat1, dat2, dat3, dat4)
    
    Sys.sleep(0.05)
    setTxtProgressBar(lbar, i)
    
  }
  
  df=df[-1,]
  #close loading bar
  close(lbar)
  print("--Data:--loading--complete--")
  
  
  lm1t1=subset(df, Mediator=="LM1T1")
  lm1t2=subset(df, Mediator=="LM1T2")
  am1t1=subset(df, Mediator=="AM1T1")
  am1t2=subset(df, Mediator=="AM1T2")
  
  
  return(list(mdf=df, lm1t1=lm1t1, lm1t2=lm1t2, am1t1=am1t1, am1t2=am1t2))
}

  
#=============================================================================================================

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
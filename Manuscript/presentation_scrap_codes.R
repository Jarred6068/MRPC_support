#============================================End--Function_1===================================================


# Lond2Addis.lookup(trio.index=613,tissue.name="WholeBlood") 

#============================================Function_2========================================================

# a function to streamline looking up trios whose classification changed from LOND ---> ADDIS
# generally used to look up individual changes after running ADDIS.PostProc
# designed to be used with the tissue names in column 2 of tissuenames.csv

Lond2Addis.lookup=function(trio.index=NULL, tissue.name=NULL, run.models=TRUE, with.pc=TRUE, correlation=TRUE, 
                           r.pack.lib="/mnt/ceph/jarredk/Rpackages"){
  
  #Load necessary files (SNP data, PC matrix and significantly associated pcs of each trio)
  snps=loadRData(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissue.name,
                       "_AllPC/data.snp.cis.trans.final.", tissue.name, ".V8.unique.snps.RData", sep = ""))
  
  pc.matrix=loadRData(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissue.name,
                            "_AllPC/PCs.matrix.",tissue.name,".RData", sep = ""))
  
  sig.pc=loadRData(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissue.name,
                         "_AllPC/List.Match.significant.trios.RData", sep = ""))
  
  library("MRPC", lib=r.pack.lib)
  
  if(with.pc==TRUE){
    
    #extract sig. asso. Pc's; get trio, and combine
    sig.asso.pcs=sig.pc[[trio.index]]
    pcs=pc.matrix[,sig.asso.pcs]
    trio=snps[,(trio.index-1)*3+(1:3)]
    trio.pc=cbind(trio, pcs)
    colnames(trio.pc)=c(colnames(trio), colnames(pc.matrix)[sig.asso.pcs])
    
    #----correlation--option----
    if(correlation==TRUE){
      r=cor(trio.pc, use = "complete.obs")
    }else{
      r=NULL
    }
    
    
    #----run--models--option----
    if(run.models==TRUE){
      
      n <- nrow (trio.pc)
      V <- colnames(trio.pc)     # Column names
      
      # Classical correlation
      suffStat <- list(C = cor(trio.pc, use = "complete.obs"),
                       n = n)
      #run MRPC on TRIO with pc's
      MRPC.fit.FDR.addis <- MRPC(trio.pc,
                                 suffStat,
                                 GV = 1,
                                 FDR = 0.05,
                                 indepTest = 'gaussCItest',
                                 labels = V,
                                 FDRcontrol = "ADDIS",
                                 verbose = FALSE)
      
      
      MRPC.fit.FDR.lond <- MRPC(trio.pc,
                                suffStat,
                                GV = 1,
                                FDR = 0.05,
                                indepTest = 'gaussCItest',
                                labels = V,
                                FDRcontrol = "LOND",
                                verbose = FALSE)
      #convert to matrix
      MRPC.fit.lond.table=as(MRPC.fit.FDR.lond@graph, "matrix")
      MRPC.fit.addis.table=as(MRPC.fit.FDR.addis@graph, "matrix")
      
      return(list(trio.with.pc=trio.pc, correlation=r, MRPC.table.LOND=MRPC.fit.lond.table, 
                  MRPC.table.ADDIS=MRPC.fit.addis.table))
      
    }else{
      
      return(list(trio.with.pc=trio.pc, correlation=r))
      
    }
    
    
    
    
    
    #--------------no---pc-----option----------------------
  }else{
    
    #get trio
    trio=snps[,(trio.index-1)*3+(1:3)]
    
    
    #----correlation--option----
    if(correlation==TRUE){
      r = cor(trio, use = "complete.obs")
    }else{
      r = NULL
    }
    
    
    #----run--models--option----
    if(run.models==TRUE){
      
      n <- nrow (trio)
      V <- colnames(trio)     # Column names
      
      # Classical correlation
      suffStat <- list(C = cor(trio, use = "complete.obs"),
                       n = n)
      #run MRPC on TRIO
      MRPC.fit.FDR.addis <- MRPC(trio,
                                 suffStat,
                                 GV = 1,
                                 FDR = 0.05,
                                 indepTest = 'gaussCItest',
                                 labels = V,
                                 FDRcontrol = "ADDIS",
                                 verbose = FALSE)
      
      
      MRPC.fit.FDR.lond <- MRPC(trio,
                                suffStat,
                                GV = 1,
                                FDR = 0.05,
                                indepTest = 'gaussCItest',
                                labels = V,
                                FDRcontrol = "LOND",
                                verbose = FALSE)
      #convert to matrix
      MRPC.fit.lond.table=as(MRPC.fit.FDR.lond@graph, "matrix")
      MRPC.fit.addis.table=as(MRPC.fit.FDR.addis@graph, "matrix")
      
      return(list(trio=trio, correlation=r, MRPC.table.LOND=MRPC.fit.lond.table, 
                  MRPC.table.ADDIS=MRPC.fit.addis.table))
      
    }else{
      
      return(list(trio=trio, correlation=r))
      
    }
    
    
  }
  
}


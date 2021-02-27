


#NOTE: tissue SmallIntestine file List.models.SmallIntestine.all.Rdata doesn't work with
#this function because it is missing data for M2



#A simple helper function for easier reading in of .Rdata files with a specific name
#================================================Helper_Function============================================

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# ADDIS.PostProc("Spleen")

#========================================START--HERE--Function_1============================================

#A function to catalog the changes in trio classification from LOND ---> ADDIS

#outputs a list of trio indices whose classifications changed
#outputs tables for each tissue showing the transitions/changes

ADDIS.PostProc=function(tissue.names=NULL){
  
  #allocate space
  tissue.list=list()
  tissue.list.diagonals=list()
  
  #create model names vec
  mod.names.Lond=c("Lond.M0","Lond.M1","Lond.M2","Lond.M3","Lond.M4")
  mod.names.Addis=c("Addis.M0","Addis.M1","Addis.M2","Addis.M3","Addis.M4")
  #loading bar
  startingbar=rep("_", length(tissue.names)*2)
  print(paste(startingbar, sep="",collapse=''))
  loading.bar=list()
  
  #looping by tissues
  for(t in 1:length(tissue.names)){
  
    #read in final output files from LOND and ADDIS versions
    file1=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/", tissue.names[t], 
                "_AllPC/List.models.", tissue.names[t], ".all.RData", sep = "")
  
    file2=paste("/mnt/ceph/jarredk/AddisReRunFiles/List.models.",
                tissue.names[t], ".all.RData", sep = "")
  
    list.Lond=loadRData(fileName=file1)
    list.Addis=loadRData(fileName=file2)  
    
    #preallocate space
    models.list=list()
    diagonals=NULL
    #first loop to find matching for M(i)
    for(j in 1:length(list.Lond)){
    
      match.trios=match(list.Lond[[j]], list.Addis[[j]])
      nomatch=list.Lond[[j]][is.na(match.trios)]
      diagonals[j]=length(match.trios[!is.na(match.trios)])
      
      location.ADDIS=list()
      #second loop to find location of missing
      for(i in 1:length(list.Lond)){
      
       trio.match.idx=match(nomatch, list.Addis[[i]])
       trio.match.idx=trio.match.idx[!is.na(trio.match.idx)]
       location.ADDIS[[i]]=list.Addis[[i]][trio.match.idx]
       
      }
    
      names(location.ADDIS)=mod.names.Addis
      models.list[[j]]=location.ADDIS
     
    
    } 
    
    names(models.list)=mod.names.Lond
    tissue.list.diagonals[[t]]=diagonals
    tissue.list[[t]]=models.list
   
    #loading bar
    loading.bar[[t]]=c(">",rep("=", t*2), startingbar[1:(length(startingbar)-t*2)])
    print(paste(paste(loading.bar[[t]], sep = "", collapse = ""), paste("...",round(((t*2)/(length(tissue.names)*2))*100,2),"%", sep = ""), sep = ""))
   
  }
  
  names(tissue.list)=tissue.names
  names(tissue.list.diagonals)=tissue.names
  
  #-------------------end of list generation-------------------
  #count up dist

  vec1=NULL
  sum2table=as.data.frame(matrix(0, nrow = 5, ncol = 5))
  all.tables=list()
  
  for(t in 1:length(tissue.names)){

    for(i in 1:length(mod.names.Lond)){

      for(j in 1:length(mod.names.Lond)){

        vec1[j]=length(tissue.list[[t]][[i]][[j]])


      }
      
      sum2table[i,]=vec1
      
    }
    
    #Add back diagonal values and naming
    diag(sum2table)=tissue.list.diagonals[[t]]
    colnames(sum2table)=mod.names.Addis
    row.names(sum2table)=mod.names.Lond
    all.tables[[t]]=sum2table
    
  }
  names(all.tables)=tissue.names

  print("done")
  
  return(output=list(tiss=tissue.list, diags=tissue.list.diagonals, all.tables=all.tables))
  
}
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
      
      return(list(trio.with.pc=trio.pc, correlation=r, MRPC.table.LOND=MRPC.fit.lond.table, 
                  MRPC.table.ADDIS=MRPC.fit.addis.table))
      
    }else{
      
      return(list(trio=trio, correlation=r))
      
    }
    
    
  }
  
}


#===================================END-Function_2=============================================


tissue.names=read.csv("/mnt/ceph/jarredk/AddisReRunFiles/tissuenames.csv", header = T)


#======================================FUCNTION_3==============================================
# this function classifies M1 trios in LOND and ADDIS as type 1 or type 2 where the definitions
# type 1: V1-->T1-->T2
# type 2: V1-->T2-->T1
# the indices for type 1 & 2 M1 trios are returned as sub-elements to the output list as $Catalog
# which can be further indexed to a specific tissue. Additionally, a table of counts and relative
# proportions is also returned

ADDIS.M1.check=function(tissue.names=NULL, verbose=FALSE){
  #SYNTAX: 
  #tissue.names -- can be either a single tissue name or character vector of names
  #corresponding to the format of column 2 in tissuenames.csv
  
  #preallocate master list
  All.Data=list()
  
  #establish ground truths for M1.1 and M1.2
  Truth.M1.type1=matrix(c(1,0,0,1,0,0), nrow = 3, ncol = 2, byrow = T)
  Truth.M1.type2=matrix(c(0,1,0,0,1,0), nrow = 3, ncol = 2, byrow = T)
  
  for(t in 1:length(tissue.names)){
    
    #preallocate list structures
    Addis.M1=list()
    Addis.M1$type1=NULL
    Addis.M1$type2=NULL
    Lond.M1=list()
    Lond.M1$type1=NULL
    Lond.M1$type2=NULL
    
    #Load needed Files
    list.Lond=loadRData(fileName=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/", tissue.names[t], 
                                       "_AllPC/List.models.", tissue.names[t], ".all.RData", sep = ""))
    
    list.Addis=loadRData(fileName=paste("/mnt/ceph/jarredk/AddisReRunFiles/List.models.",
                                        tissue.names[t], ".all.RData", sep = ""))
    #get M1's from ADDIS/LOND
    M1s.addis=list.Addis$M1
    M1s.lond=list.Lond$M1
    #loop for LOND
    for(i in 1:length(M1s.lond)){
      
      results.lond=Lond2Addis.lookup(trio.index = M1s.lond[i], tissue.name = tissue.names[t])
      MRPC.lond=results.lond$MRPC.table.LOND[1:3,2:3]
      #print(MRPC.lond)
      
      Lond.M1$type1[i]=ifelse(all.equal(MRPC.lond, Truth.M1.type1, check.attributes=FALSE)==TRUE, M1s.lond[i], NA)
      #print(Lond.M1$type1)
      Lond.M1$type2[i]=ifelse(all.equal(MRPC.lond, Truth.M1.type2, check.attributes=FALSE)==TRUE, M1s.lond[i], NA)
      #print(Lond.M1$type2)
      
    }
    
    
    #loop for ADDIS
    for(i in 1:length(M1s.addis)){
      
      results.addis=Lond2Addis.lookup(trio.index = M1s.addis[i], tissue.name = tissue.names[t])
      MRPC.addis=results.addis$MRPC.table.ADDIS[1:3,2:3]
      
      Addis.M1$type1[i]=ifelse(all.equal(MRPC.addis, Truth.M1.type1, check.attributes=FALSE)==TRUE, M1s.addis[i], NA)
      Addis.M1$type2[i]=ifelse(all.equal(MRPC.addis, Truth.M1.type2, check.attributes=FALSE)==TRUE, M1s.addis[i], NA)
      
    }
    
    #Printing to make sure M1.1+M1.2 = M1 total
    if(verbose==TRUE){
      print(c(sum(length(Addis.M1$type1[!is.na(Addis.M1$type1)]), length(Addis.M1$type2[!is.na(Addis.M1$type2)])), length(M1s.addis)))
      print(c(sum(length(Lond.M1$type1[!is.na(Lond.M1$type1)]), length(Lond.M1$type2[!is.na(Lond.M1$type2)])), length(M1s.lond)))
    }
    
    #Remove NA's
    Lond.M1$type1=Lond.M1$type1[!is.na(Lond.M1$type1)]
    Lond.M1$type2=Lond.M1$type2[!is.na(Lond.M1$type2)]
    Addis.M1$type1=Addis.M1$type1[!is.na(Addis.M1$type1)]
    Addis.M1$type2=Addis.M1$type2[!is.na(Addis.M1$type2)]
    
    #Add to master list
    All.Data$tiss[[t]]=list(LOND=Lond.M1, ADDIS=Addis.M1)
    
  }
  
  names(All.Data$tiss)=tissue.names
  #preallocate
  vec1=NULL
  vec2=NULL
  df.count=as.data.frame(matrix(0, nrow = length(All.Data$tiss), ncol = 6))
  df.per=as.data.frame(matrix(0, nrow = length(All.Data$tiss), ncol = 6))
  #calc summary percentages and totals
  for(j in 1:length(All.Data$tiss)){
    
    list.Lond=loadRData(fileName=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/", tissue.names[j], 
                                       "_AllPC/List.models.", tissue.names[j], ".all.RData", sep = ""))
    
    list.Addis=loadRData(fileName=paste("/mnt/ceph/jarredk/AddisReRunFiles/List.models.",
                                        tissue.names[j], ".all.RData", sep = ""))
    M1s.addis=list.Addis$M1
    M1s.lond=list.Lond$M1
    
    total.lond=length(M1s.lond)
    total.addis=length(M1s.addis)
    
    counts.addis=c(length(All.Data$tiss[[j]]$ADDIS$type1), length(All.Data$tiss[[j]]$ADDIS$type2))
    counts.lond=c(length(All.Data$tiss[[j]]$LOND$type1), length(All.Data$tiss[[j]]$LOND$type2))
    
    vec1=c(counts.addis, total.addis, counts.lond,total.lond)
    vec2=c(counts.addis/total.addis, total.addis, counts.lond/total.lond, total.lond)
    
    df.count[j, ] = vec1
    df.per[j, ] = vec2
    
  }
  
  #naming
  colnames(df.count)=c("ADDIS.M1.1","ADDIS.M1.2","ADDIS.TOTAL", "LOND.M1.1","LOND.M1.2", "LOND.TOTAL")
  row.names(df.count)=tissue.names
  colnames(df.per)=c("%_ADDIS.M1.1","%_ADDIS.M1.2","ADDIS.TOTAL", "%_LOND.M1.1","%_LOND.M1.2","LOND.TOTAL")
  row.names(df.per)=tissue.names
  
  
  return(list(Catalog=All.Data, counts=df.count, relative.percentages=df.per))
  
}

#==========================================END--FUNCTION_3===================================
  
  
  
  
  
  
  
  
  














































#read in the trio specific algorithm functions
source("/mnt/ceph/jarredk/MRPC_UPDATE/mrnet_triospecific.R")

##########################################################
#this function applies the mrnet trio specific algorithm to
#a single tissue specied by the tissue name 

mrnet.analyze=function(GTEx.tissname="WholeBlood", save=TRUE, verbose=TRUE){
  
  

  #load in useful functions from ADDIS_verify
  source("/mnt/ceph/jarredk/ADDIS_verify/ADDIS_Post_Analysis_processing.R")
  #Load necessary files (trio data, PC matrix and significantly associated pcs of each trio)
  trios=loadRData(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",GTEx.tissname,
                        "_AllPC/data.snp.cis.trans.final.", GTEx.tissname, ".V8.unique.snps.RData", sep = ""))
  
  pc.matrix=loadRData(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",GTEx.tissname,
                            "_AllPC/PCs.matrix.",GTEx.tissname,".RData", sep = ""))
  
  sig.pc=loadRData(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",GTEx.tissname,
                         "_AllPC/List.Match.significant.trios.RData", sep = ""))
  
  #get the number of trios for tissue + indices
  num.of.trios=dim(trios)[2]/3
  trio.index=c(1:num.of.trios)
  trio.list=list()
  
  #allocate all trios to a list
  for(i in 1:num.of.trios){
    #extract sig. asso. Pc's; get trio, and combine
    sig.asso.pcs=sig.pc[[trio.index[i]]]
    pcs=pc.matrix[,sig.asso.pcs]
    trio.with.pc=cbind.data.frame(trios[,(trio.index[i]-1)*3+(1:3)],pcs)
    if(length(sig.asso.pcs)==0){
      colnames(trio.with.pc)=c(colnames(trios[,(trio.index[i]-1)*3+(1:3)]))
    }else{
      colnames(trio.with.pc)=c(colnames(trios[,(trio.index[i]-1)*3+(1:3)]),paste0("PC",sig.asso.pcs))
    }
    
    trio.list[[i]]=trio.with.pc
    
  }
  
  #preform regressions and classify model types
  reg.res=list()
  inf.mods=list()
  
    
  reg.res=sapply(trio.list, infer.trio)
  inf.mods=apply(reg.res, 2, class.vec)
  
  final=summary(as.factor(inf.mods))
  
  if(save==TRUE){
    save(reg.res, file = "/mnt/ceph/jarredk/MRPC_UPDATE/GTEx_tissues_analysis/WholeBlood/mrnet.infdata.WB.RData")
    save(inf.mods, file = "/mnt/ceph/jarredk/MRPC_UPDATE/GTEx_tissues_analysis/WholeBlood/mrnet.inf.mods.WB.RData")
    save(final, file = "/mnt/ceph/jarredk/MRPC_UPDATE/GTEx_tissues_analysis/WholeBlood/mrnet.final.RData")
  }
  
  return(final)
}






























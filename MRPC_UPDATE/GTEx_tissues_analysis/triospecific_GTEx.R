
#read in the trio specific algorithm functions
source("/mnt/ceph/jarredk/MRPC_UPDATE/mrnet_triospecific.R")
#load in useful functions from ADDIS_verify
source("/mnt/ceph/jarredk/ADDIS_verify/ADDIS_Post_Analysis_processing.R")

##########################################################
#this function applies the mrnet trio specific algorithm to
#a single tissue specied by the tissue name 

mrnet.analyze=function(GTEx.tissname="WholeBlood", save=TRUE, verbose=TRUE){
  
  final=vector("list", length = length(GTEx.tissname))
  names(final)=GTEx.tissname
  
  for(t in 1:length(GTEx.tissname)){
    #Load necessary files (trio data, PC matrix and significantly associated pcs of each trio)
    trios=loadRData(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",GTEx.tissname[t],
                          "_AllPC/data.snp.cis.trans.final.", GTEx.tissname[t], ".V8.unique.snps.RData", sep = ""))
    
    pc.matrix=loadRData(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",GTEx.tissname[t],
                              "_AllPC/PCs.matrix.",GTEx.tissname[t],".RData", sep = ""))
    
    sig.pc=loadRData(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",GTEx.tissname[t],
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
    
    #trio.list=trio.list[1:5]
    #print(lapply(trio.list, head))
    #preform regressions and classify model types
    reg.res=list()
    inf.mods=list()
    
    
    reg.res=sapply(trio.list, infer.trio)
    #print(reg.res)
    inf.mods=apply(reg.res, 2, class.vec)
    #print(inf.mods)
    
    sum.inf=summary(as.factor(inf.mods))
    #print(sum.inf)
    
    if(save==TRUE){
      save(reg.res, file = paste0("/mnt/ceph/jarredk/MRPC_UPDATE/GTEx_tissues_analysis/", GTEx.tissname[t], "/mrnet.infdata.WB.RData"))
      save(inf.mods, file = paste0("/mnt/ceph/jarredk/MRPC_UPDATE/GTEx_tissues_analysis/", GTEx.tissname[t],"/mrnet.inf.mods.WB.RData"))
      save(sum.inf, file = paste0("/mnt/ceph/jarredk/MRPC_UPDATE/GTEx_tissues_analysis/", GTEx.tissname[t],"/mrnet.final.RData"))
    }
    
    final[[t]]=sum.inf
    
  }
  
  return(final)
}






























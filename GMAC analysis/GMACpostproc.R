
source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")

library('qvalue', lib="/mnt/ceph/jarredk/Rpackages")
library('MRPC', lib="/mnt/ceph/jarredk/Rpackages")

source("/mnt/ceph/jarredk/ADDIS_verify/ADDIS_Post_Analysis_processing.R")

top5=c(1,6, 33, 40, 48)
path='/mnt/ceph/jarredk/Reg_Net/'
tissues.vec=tissue.names[top5, 2:3]




#function to run post processing on GMAC output

run.postproc=function(output.pvals=NULL, trio.ref=NULL){
  
  qvals=qvalue(output.pvals, fdr.level = 0.1)
  
  num.inferred=sum(qvals$significant)
  
  sig.trios=trio.ref[qvals$significant,]
  
  result=cbind.data.frame(qvals$qvalue[qvals$significant], output.pvals[qvals$significant])
  
  
  return(list(total.inferred=num.inferred, sig.trios=sig.trios, pq.values=result))
}

# #load in output for the cis trios:
# 
# output.WB=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/WholeBlood/all_trios_output_cis.Rdata')
# output.AS=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/AdiposeSubcutaneous/all_trios_output_cis.Rdata')
# output.ArT=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/ArteryTibial/all_trios_output_cis.Rdata')
# output.MS=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/MuscleSkeletal/all_trios_output_cis.Rdata')
# output.SSE=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/SkinSunExposed/all_trios_output_cis.Rdata')
# 
# #post process GMAC-cis-mediated
# result.WB=run.postproc(output.pvals = output.WB$output.table[,5], trio.ref = output.WB$output.table[,1:3])
# result.AS=run.postproc(output.pvals = output.AS$output.table[,5], trio.ref = output.AS$output.table[,1:3])
# result.ArT=run.postproc(output.pvals = output.ArT$output.table[,5], trio.ref = output.ArT$output.table[,1:3])
# result.MS=run.postproc(output.pvals = output.MS$output.table[,5], trio.ref = output.MS$output.table[,1:3])
# result.SSE=run.postproc(output.pvals = output.SSE$output.table[,5], trio.ref = output.SSE$output.table[,1:3])
# 
# 
# #read in output for the trans trios:
# 
# #load in output for the cis trios:
# 
# output.t.WB=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/WholeBlood/all_trios_output_trans.Rdata')
# output.t.AS=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/AdiposeSubcutaneous/all_trios_output_trans.Rdata')
# output.t.ArT=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/ArteryTibial/all_trios_output_trans.Rdata')
# output.t.MS=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/MuscleSkeletal/all_trios_output_trans.Rdata')
# output.t.SSE=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/SkinSunExposed/all_trios_output_trans.Rdata')
# 
# result.t.WB=run.postproc(output.pvals = output.t.WB$output.table[,5], trio.ref = output.t.WB$output.table[,1:3])
# result.t.AS=run.postproc(output.pvals = output.t.AS$output.table[,5], trio.ref = output.t.AS$output.table[,1:3])
# result.t.ArT=run.postproc(output.pvals = output.t.ArT$output.table[,5], trio.ref = output.t.ArT$output.table[,1:3])
# result.t.MS=run.postproc(output.pvals = output.t.MS$output.table[,5], trio.ref = output.t.MS$output.table[,1:3])
# result.t.SSE=run.postproc(output.pvals = output.t.SSE$output.table[,5], trio.ref = output.t.SSE$output.table[,1:3])
# 
# 
# result.table=cbind.data.frame(tissues.vec[,1],
#                               c(dim(output.AS$output.table)[1],
#                                 dim(output.ArT$output.table)[1],
#                                 dim(output.MS$output.table)[1],
#                                 dim(output.SSE$output.table)[1],
#                                 dim(output.WB$output.table)[1]),
#                               c(result.AS$total.inferred,
#                                 result.ArT$total.inferred,
#                                 result.MS$total.inferred,
#                                 result.SSE$total.inferred,
#                                 result.WB$total.inferred),
#                               c(result.t.AS$total.inferred,
#                                 result.t.ArT$total.inferred,
#                                 result.t.MS$total.inferred,
#                                 result.t.SSE$total.inferred,
#                                 result.t.WB$total.inferred),
#                               c(123, 90, 104, 122, 131),
#                               c(130, 113, 103, 123, 113))
# 
# 
# colnames(result.table)=c("Tissue",
#                          "Total.Num.Samples",
#                          "GMAC.Total.M1T1",
#                          "GMAC.Total.M1T2",
#                          "LOND.total.M1",
#                          "ADDIS.total.M1")



write.csv(result.table, file="/mnt/ceph/jarredk/GMACanalysis/final_result_table.csv")



match.trios=function(tissues=tissues.vec[,1], which.mrpc="ADDIS"){
  
  library('prodlim', lib="/mnt/ceph/jarredk/Rpackages")
  matched.cis=NULL
  matched.trans=NULL
  
  unmatched.cis=vector("list", length = length(tissues))
  unmatched.trans=vector("list", length = length(tissues))
  
  for(i in 1:length(tissues)){
    
    
    #read in GMAC results
    output.cis=loadRData(fileName=paste('/mnt/ceph/jarredk/GMACanalysis/', tissues[i], '/all_trios_output_cis50.Rdata', sep = ""))
    output.trans=loadRData(fileName=paste('/mnt/ceph/jarredk/GMACanalysis/', tissues[i], '/all_trios_output_trans50.Rdata', sep=""))
    
    #process 
    postproc.cis=run.postproc(output.pvals = output.cis$output.table[,5], 
                            trio.ref = output.cis$output.table[,1:3])
    
    postproc.trans=run.postproc(output.pvals = output.trans$output.table[,5], 
                              trio.ref = output.trans$output.table[,1:3])
    
    if(which.mrpc=="ADDIS"){
      
      m1t1=read.csv(file = paste("/mnt/ceph/jarredk/Reg_Net/AM1T1/", tissues[i], ".csv", sep = ""))
      m1t2=read.csv(file = paste("/mnt/ceph/jarredk/Reg_Net/AM1T2/", tissues[i], ".csv", sep = ""))
      
    }else{
      
      m1t1=read.csv(file = paste("/mnt/ceph/jarredk/Reg_Net/LM1T1/", tissues[i], ".csv", sep = ""))
      m1t2=read.csv(file = paste("/mnt/ceph/jarredk/Reg_Net/LM1T2/", tissues[i], ".csv", sep = ""))
      
    }
    
    
    matched.cis[i]=length(na.omit(row.match(m1t1[,2:4],postproc.cis$sig.trios)))
    matched.trans[i]=length(na.omit(row.match(m1t2[,2:4],postproc.trans$sig.trios)))
    
    unmatched.cis[[i]]=postproc.cis$sig.trios[-c(na.omit(row.match(m1t1[,2:4],postproc.cis$sig.trios))),]
    unmatched.trans[[i]]=postproc.trans$sig.trios[-c(na.omit(row.match(m1t2[,2:4],postproc.trans$sig.trios))),]
    

    
  }
  
  names(unmatched.cis)=tissues
  names(unmatched.trans)=tissues
  
  table1=cbind.data.frame(tissues, matched.cis, matched.trans)
  colnames(table1)=c("tissue","cis.common", "trans.common")
  
  return(list(table=table1, unmatched.cis=unmatched.cis, unmatched.trans=unmatched.trans))
    
    
}












#match.trios(tissues=tissues.vec[,1], which.mrpc="ADDIS")
#match.trios(tissues=tissues.vec[,1], which.mrpc="LOND")
  


#a simple function to find the trio index of a trio based on the SNP name

locate.trio=function(input=NULL, trio.mat=NULL){
  
  index.vec=NULL
  
  for(j in 1:dim(input)[1]){
    
    ind1=which(colnames(trio.mat)==input[j,1])
    index.vec[j]=(ind1-1)/3+1
    
  }
  
  return(index.vec)
  
}

# a simple logical function passed to lapply in which.class()
fn=function(x,b){x==b}

#a simple function to find the addis classification of a trio

which.class=function(ind=NULL, classes.list=NULL){
  
  models=c("MO","M1","M2","M3","M4")
  class.vec=NULL
  
  for(i in 1:length(ind)){
    
    vec1=unlist(lapply(lapply(classes.list, fn, ind[i]),any))
    
    class.vec[i]=ifelse(any(vec1), models[vec1], "Other")
    
  }
  
  return(class.vec)
  
}



#cross analyze mediation test results with ADDIS trios
cross.analyze=function(tissues=tissues.vec[,1], which.type="cis", save=FALSE, 
                       path="/mnt/ceph/jarredk/GMACanalysis/"){
  
  out1=match.trios(tissues=tissues, which.mrpc="ADDIS")
  final.list=vector("list", length = length(tissues))
  
  for(i in 1:length(tissues)){
    
    addis.classes=loadRData(fileName=paste("/mnt/ceph/jarredk/AddisReRunFiles/List.models.", tissues[i], ".all.RData" ,sep=""))
    
      
    file1=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",
                tissues[i],"_AllPC/data.snp.cis.trans.final.",
                tissues[i],".V8.unique.snps.RData", sep = "")
    
    trios=loadRData(fileName=file1)
    
    if(which.type=="cis"){
      
      indices=locate.trio(out1$unmatched.cis[[i]], trios)
      classes=which.class(ind = indices, classes.list = addis.classes)
      tab=cbind.data.frame(out1$unmatched.cis[[i]], indices, classes)
      colnames(tab)=c(colnames(out1$unmatched.cis[[i]]), "index.in.trio.mat","Addis.infer.Class")
      final.list[[i]]=tab
      
    }else{
      
      indices=locate.trio(out1$unmatched.trans[[i]], trios)
      classes=which.class(ind = indices, classes.list = addis.classes)
      tab=cbind.data.frame(out1$unmatched.trans[[i]], indices, classes)
      colnames(tab)=c(colnames(out1$unmatched.trans[[i]]), "index.in.trio.mat","Addis.infer.Class")
      final.list[[i]]=tab
      
    }
    
    
    
    
    
  }
  
  if(save==TRUE){
    
    save(final.list, file = paste0(path, "GMAC.to.Addis.Class.List.",which.type, ".RData"))
    
  }
  
  return(final.list)
  
}


#plots the distribution of M types for non-overlapping trios classified as mediation in GMAC

plot.dist=function(class.vec.cis=NULL, class.vec.trans=NULL){
  
  D1=summary(as.factor(class.vec.cis))
  print(D1)
  D1.1=cbind.data.frame(names(D1),D1, rep("cis", length(D1)))
  colnames(D1.1)=c("Class", "Number of Trios", "mediator type")
  
  D2=summary(as.factor(class.vec.trans))
  print(D2)
  D2.1=cbind.data.frame(names(D2),D2, rep("trans", length(D2)))
  colnames(D2.1)=c("Class", "Number of Trios", "mediator type")
  
  all.data=rbind.data.frame(D1.1,D2.1)
  all.data$Class=as.factor(all.data$Class)
  all.data$`mediator type`=as.factor(all.data$`mediator type`)
  
  library(ggpubr)
  
  
  ggbarplot(all.data, x = "Class", y = "Number of Trios",
            fill = "mediator type", color = "mediator type", 
            palette = c("cadetblue1","red4", "yellow","navyblue","salmon1", "chartreuse3"),
            label = FALSE, ylab = FALSE, position=position_dodge(0.9))+
    theme(axis.text.x = element_text(size=10, hjust=0.5,vjust=0.7), 
          axis.text.y = element_text(size=9))+
    scale_y_continuous(breaks = seq(0, 500, 50), lim = c(0, 500))+
    coord_flip()+
    ggtitle("ADDIS Model Classification for Non-Overlapping M1 Trios", 
            subtitle = "Inferred by GMAC analysis")
  
  
}







#a function to look at the regression on the trans gene in GMAC


cross.regress=function(tissue="WholeBlood", trio.ind=NULL, mod.type="cis", addis.pcs=NULL){
  
  if(mod.type=="cis"){
    out.data=loadRData(fileName=paste0('/mnt/ceph/jarredk/GMACanalysis/', tissue, '/all_trios_output_cis.Rdata'))
    
  }else{
    
    out.data=loadRData(fileName=paste0('/mnt/ceph/jarredk/GMACanalysis/', tissue, '/all_trios_output_trans.Rdata'))
    
  }
  
  trio.loc=as.matrix(out.data$input.list$trios.idx)[trio.ind,]
  SNP=t(as.data.frame(out.data$input.list$snp.dat.cis)[trio.loc[1],])
  GE=t(as.data.frame(out.data$input.list$exp.dat)[trio.loc[-1],])
  which.covs=which(out.data$cov.indicator.list[trio.ind,]==1)
  sel.cov.pool=as.data.frame(out.data$input.list$cov.pool)[which.covs,]
  
  print("-----------GMAC-Selected-PCs--------------")
  print(row.names(sel.cov.pool))
  
  sel.cov.pool=t(sel.cov.pool)
  
  data.mat=cbind.data.frame(GE, SNP, sel.cov.pool)
  which.covs.addis=match(addis.pcs, row.names(out.data$input.list$cov.pool))
  addis.data=cbind.data.frame(GE, SNP, t(out.data$input.list$cov.pool[which.covs.addis,]))
  
  if(mod.type=="cis"){
    
    colnames(data.mat)=c("cis.gene", "trans.gene", "SNP", colnames(sel.cov.pool))
    colnames(addis.data)=c("cis.gene", "trans.gene", "SNP", addis.pcs)
    
  }else{
    
    colnames(data.mat)=c("trans.gene", "cis.gene", "SNP", colnames(sel.cov.pool))
    colnames(addis.data)=c("trans.gene", "cis.gene", "SNP", addis.pcs)
    
  }
  

  
  
  
  model.GMAC=lm(trans.gene~., data = data.mat)
  print("--------------------------GMAC------------------------------")
  print(summary(model.GMAC))
  model.ADDIS=lm(trans.gene~., data = addis.data)
  print("--------------------------ADDIS------------------------------")
  print(summary(model.ADDIS))
  
  
}






































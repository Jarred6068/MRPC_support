
source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")

library('qvalue', lib="/mnt/ceph/jarredk/Rpackages")
library('MRPC', lib="/mnt/ceph/jarredk/Rpackages")
library('permute', lib="/mnt/ceph/jarredk/Rpackages")

source("/mnt/ceph/jarredk/ADDIS_verify/ADDIS_Post_Analysis_processing.R")

top5=c(1, 6, 33, 40, 48)
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
output.WB=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/WholeBlood/all_trios_output_cis.Rdata')
output.AS=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/AdiposeSubcutaneous/all_trios_output_cis.Rdata')
output.ArT=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/ArteryTibial/all_trios_output_cis.Rdata')
output.MS=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/MuscleSkeletal/all_trios_output_cis.Rdata')
output.SSE=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/SkinSunExposed/all_trios_output_cis.Rdata')
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
output.t.WB=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/WholeBlood/all_trios_output_trans.Rdata')
output.t.AS=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/AdiposeSubcutaneous/all_trios_output_trans.Rdata')
output.t.ArT=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/ArteryTibial/all_trios_output_trans.Rdata')
output.t.MS=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/MuscleSkeletal/all_trios_output_trans.Rdata')
output.t.SSE=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/SkinSunExposed/all_trios_output_trans.Rdata')
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



#write.csv(result.table, file="/mnt/ceph/jarredk/GMACanalysis/final_result_table.csv")



#a function to match which trios are found in the AM1T1/T2 and LM1T1/T2 tables and which have 
#inferred model types that differ

match.trios=function(tissues=tissues.vec[,1], which.mrpc="ADDIS"){
  
  library('prodlim', lib="/mnt/ceph/jarredk/Rpackages")
  
  num.matched.cis=NULL
  num.matched.trans=NULL
  
  matched.cis=vector("list", length = length(tissues))
  matched.trans=vector("list", length = length(tissues))
  
  unmatched.cis=vector("list", length = length(tissues))
  unmatched.trans=vector("list", length = length(tissues))
  
  for(i in 1:length(tissues)){
    
    
    #read in GMAC results
    output.cis=loadRData(fileName=paste('/mnt/ceph/jarredk/GMACanalysis/', tissues[i], '/all_trios_output_cis.Rdata', sep = ""))
    output.trans=loadRData(fileName=paste('/mnt/ceph/jarredk/GMACanalysis/', tissues[i], '/all_trios_output_trans.Rdata', sep=""))
    
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
    
    
    num.matched.cis[i]=length(na.omit(row.match(m1t1[,2:4],postproc.cis$sig.trios)))
    num.matched.trans[i]=length(na.omit(row.match(m1t2[,2:4],postproc.trans$sig.trios)))
    
    unmatched.cis[[i]]=postproc.cis$sig.trios[-c(na.omit(row.match(m1t1[,2:4],postproc.cis$sig.trios))),]
    unmatched.trans[[i]]=postproc.trans$sig.trios[-c(na.omit(row.match(m1t2[,2:4],postproc.trans$sig.trios))),]
    
    matched.cis[[i]]=postproc.cis$sig.trios[c(na.omit(row.match(m1t1[,2:4],postproc.cis$sig.trios))),]
    matched.trans[[i]]=postproc.trans$sig.trios[c(na.omit(row.match(m1t2[,2:4],postproc.trans$sig.trios))),]
    

    
  }
  
  names(unmatched.cis)=tissues
  names(unmatched.trans)=tissues
  
  names(matched.cis)=tissues
  names(matched.trans)=tissues
  
  table1=cbind.data.frame(tissues, num.matched.cis, num.matched.trans)
  colnames(table1)=c("tissue","cis.common", "trans.common")
  
  return(list(table=table1, 
              unmatched.cis=unmatched.cis, 
              unmatched.trans=unmatched.trans,
              matched.cis=matched.cis,
              matched.trans=matched.trans))
    
    
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
  
  
  final.list.unmatched=vector("list", length = length(tissues))
  final.list.matched=vector("list", length = length(tissues))
  
  for(i in 1:length(tissues)){
    
    out1=match.trios(tissues=tissues[i], which.mrpc="ADDIS")
    
    addis.classes=loadRData(fileName=paste("/mnt/ceph/jarredk/AddisReRunFiles/List.models.", tissues[i], ".all.RData" ,sep=""))
    
      
    file1=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",
                tissues[i],"_AllPC/data.snp.cis.trans.final.",
                tissues[i],".V8.unique.snps.RData", sep = "")
    
    trios=loadRData(fileName=file1)
    
    if(which.type=="cis"){
      
      indices=locate.trio(out1$unmatched.cis[[1]], trios)
      indices2=locate.trio(out1$matched.cis[[1]], trios)
      classes=which.class(ind = indices, classes.list = addis.classes)
      classes2=which.class(ind = indices2, classes.list = addis.classes)
      tab=cbind.data.frame(out1$unmatched.cis[[1]], indices, classes)
      tab2=cbind.data.frame(out1$matched.cis[[1]], indices2, classes2)
      colnames(tab)=c(colnames(out1$unmatched.cis[[1]]), "index.in.trio.mat","Addis.infer.Class")
      colnames(tab2)=c(colnames(out1$matched.cis[[1]]), "index.in.trio.mat","Addis.infer.Class")
      final.list.unmatched[[i]]=tab
      final.list.matched[[i]]=tab2
      
    }else{
      
      indices=locate.trio(out1$unmatched.trans[[i]], trios)
      indices2=locate.trio(out1$matched.trans[[i]], trios)
      classes=which.class(ind = indices, classes.list = addis.classes)
      classes2=which.class(ind = indices2, classes.list = addis.classes)
      tab=cbind.data.frame(out1$unmatched.trans[[i]], indices, classes)
      tab2=cbind.data.frame(out1$matched.trans[[i]], indices2, classes2)
      colnames(tab)=c(colnames(out1$unmatched.trans[[i]]), "index.in.trio.mat","Addis.infer.Class")
      colnames(tab2)=c(colnames(out1$matched.trans[[i]]), "index.in.trio.mat","Addis.infer.Class")
      final.list.unmatched[[i]]=tab
      final.list.matched[[i]]=tab2
      
    }
    
    
    
    
    
  }
  
  if(save==TRUE){
    
    save(final.list.unmatched, file = paste0(path, "GMAC.to.Addis.Class.List.Unmatched.",which.type, ".RData"))
    save(final.list.matched, file = paste0(path, "GMAC.to.Addis.Class.List.Matched",which.type, ".RData"))
    
  }
  
  return(list(fl.unm=final.list.unmatched, fl.m=final.list.matched))
  
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
  known.conf=t(out.data$input.list$known.conf)
  SNP=t(as.data.frame(out.data$input.list$snp.dat.cis)[trio.loc[1],])
  GE=t(as.data.frame(out.data$input.list$exp.dat)[trio.loc[-1],])
  which.covs=which(out.data$cov.indicator.list[trio.ind,]==1)
  sel.cov.pool=as.data.frame(out.data$input.list$cov.pool)[which.covs,]
  
  print("-----------GMAC-Selected-PCs--------------")
  print(row.names(sel.cov.pool))
  
  sel.cov.pool=t(sel.cov.pool)
  
  data.mat=cbind.data.frame(GE, SNP, sel.cov.pool, known.conf)
  which.covs.addis=match(addis.pcs, row.names(out.data$input.list$cov.pool))
  addis.data=cbind.data.frame(GE, SNP, t(out.data$input.list$cov.pool[which.covs.addis,]))
  
  if(mod.type=="cis"){
    
    colnames(data.mat)=c("cis.gene", "trans.gene", "SNP", colnames(sel.cov.pool), colnames(known.conf))
    colnames(addis.data)=c("cis.gene", "trans.gene", "SNP", addis.pcs)
    
  }else{
    
    colnames(data.mat)=c("trans.gene", "cis.gene", "SNP", colnames(sel.cov.pool), colnames(known.conf))
    colnames(addis.data)=c("trans.gene", "cis.gene", "SNP", addis.pcs)
    
  }
  
  data.mat2=data.mat
  addis.data2=addis.data
  #convert.factors
  #data.mat$SNP=as.factor(data.mat$SNP)
  data.mat$pcr=as.factor(data.mat$pcr)
  data.mat$sex=as.factor(data.mat$sex)
  data.mat$platform=as.factor(data.mat$platform)
  #addis.data$SNP=as.factor(addis.data$SNP)
  
  model.GMAC=lm(trans.gene~., data = data.mat)
  print("--------------------------GMAC------------------------------")
  print(summary(model.GMAC))
  model.ADDIS=lm(trans.gene~., data = addis.data)
  print("--------------------------ADDIS------------------------------")
  print(summary(model.ADDIS))
  
  model.matrix.addis=model.matrix(model.ADDIS)
  model.matrix.gmac=model.matrix(model.GMAC)
  b.gmac=summary(model.GMAC)$coefficients[,1]
  b.addis=summary(model.ADDIS)$coefficients[,1]
  sigma.gmac=summary(model.GMAC)$sigma
  sigma.addis=summary(model.ADDIS)$sigma
  
  return(list(addis=addis.data2, 
              GMAC=data.mat2, 
              Regress=list(X.addis=model.matrix.addis,
                           X.gmac=model.matrix.gmac,
                           b.gmac=b.gmac,
                           b.addis=b.addis,
                           sigma.gmac=sigma.gmac,
                           sigma.addis=sigma.addis)))
  

  
}


#a function to recreate the permutation regression from GMAC

run.permuted.reg=function(trio, nperms=1000, plot=TRUE, filename=NULL){
  
  wald.stat=NULL
  
  trio2=trio
  #convert.factors
  #trio2$SNP=as.factor(trio2$SNP)
  trio2$pcr=as.factor(trio2$pcr)
  trio2$sex=as.factor(trio2$sex)
  trio2$platform=as.factor(trio2$platform)
  
  test.wald=summary(lm(trans.gene~., data=trio2))$coefficients[2,3]
  
  for(i in 1:nperms){
    
    #conditioning on each genotype
    trio0=trio[which(trio$SNP==0),]
    trio1=trio[which(trio$SNP==1),]
    trio2=trio[which(trio$SNP==2),]
    
    #permute
    trio0$cis.gene=trio0$cis.gene[shuffle(trio0$cis.gene)]
    trio1$cis.gene=trio1$cis.gene[shuffle(trio1$cis.gene)]
    trio2$cis.gene=trio2$cis.gene[shuffle(trio2$cis.gene)]
    
    #print(trio0$cis.gene[shuffle(trio0$cis.gene)][1:5])
    
    trio.permuted=rbind.data.frame(trio0, trio1, trio2)
    #print(i)
    #print(head(trio.permuted))
    
    #trio.permuted$SNP=as.factor(trio.permuted$SNP)
    trio.permuted$pcr=as.factor(trio.permuted$pcr)
    trio.permuted$sex=as.factor(trio.permuted$sex)
    trio.permuted$platform=as.factor(trio.permuted$platform)
    
    wald.stat[i]=summary(lm(trans.gene~., data=trio.permuted))$coefficients[2,3]
    
    
  }
  
  cmap=c(rep(0, dim(trio0)[1]), rep(1, dim(trio1)[1]), rep(2, dim(trio2)[1]))
  
  if(plot==TRUE){
    
    
    png(filename = filename)
    
    par(mfrow=c(2,2))
    plot(trio.permuted$cis.gene,
         trio.permuted$trans.gene,
         xlab="Permuted Cis Gene",
         ylab="Trans Gene",
         col=as.factor(cmap), 
         bg=as.factor(cmap), 
         pch=21)
    
    legend("topleft", c("Gt=0", "Gt=1","Gt=2"), col = c(1,2,3), pch=c(21,21,21), fill=c(1,2,3))
    
    plot(trio$cis.gene,
         trio$trans.gene,
         xlab="Cis Gene",
         ylab="Trans Gene",
         col=as.factor(trio$SNP), 
         bg=as.factor(trio$SNP), 
         pch=21)
    
    hist(wald.stat)
    
    par(mfrow=c(1,1))
    
    dev.off()
    
    
    
  }
  
  p.value=sum(abs(wald.stat)>abs(test.wald))/length(wald.stat)
  
  return(list(p.value=p.value, null.wald.stats=wald.stat, obs.wald.stat=test.wald))
  
}



#simple function used in check.ns.gmac() for coefficient
#simulation matrix to longform for plotting

longform=function(iters=NULL,varvec1=NULL, cutdata=NULL){
  
  var1=rep(varvec1, iters)
  
  long=cbind.data.frame(Coef=var1, Value=as.vector(as.matrix(cutdata)))
  return(long)
  
}

#a function to check the numerical stability of the regression using all adaptively selected pcs
#in GMAC, and using the coefficients and residual standard error obtained from the regression of a 
#specific trio (the output of cross.regress). The trans gene is simulated by randomly sampling the noise
#at each iteration and running the regression on the simulated trans gene

#output: indicator matrix for which regressions had significant SNP (col 1) and significant cis.gene (col 2)

check.ns.gmac=function(input.list.data=NULL, iters=1, alpha=0.01, print.it=FALSE, which.mod="GMAC",
                       plot.it=FALSE, save.name=NULL){
  
  sig.snp=NULL
  sig.cis=NULL
  if(which.mod=="GMAC"){
    
    X=as.matrix(input.list.data$Regress$X.gmac)
    b=input.list.data$Regress$b.gmac
    sigma=input.list.data$Regress$sigma.gmac
    
  }else{
    
    X=as.matrix(input.list.data$Regress$X.addis)
    b=input.list.data$Regress$b.addis
    sigma=input.list.data$Regress$sigma.addis
    
  }

  simu.coefs=as.data.frame(matrix(0, nrow=length(b), ncol=iters))
  row.names(simu.coefs)=names(b)
  colnames(simu.coefs)=paste0("iter", c(1:iters))
  
  for(i in 1:iters){
    
    if(which.mod=="GMAC"){
      new.data=input.list.data$GMAC
    }else{
      new.data=input.list.data$addis
    }
    
    trans.gene.new=X%*%b+rnorm(dim(X)[1], mean = 0, sd = sigma)
    #print(trans.gene.new)
    
    new.data$trans.gene=trans.gene.new
    if(which.mod=="GMAC"){
      new.data$pcr=as.factor(new.data$pcr)
      new.data$sex=as.factor(new.data$sex)
      new.data$platform=as.factor(new.data$platform)
      #new.data$SNP=as.factor(new.data$SNP)
    }
    
    reg=summary(lm(trans.gene~., data = new.data))
    if(print.it==TRUE){print(reg)}
    coef=as.data.frame(reg$coefficient)
    p.snp=coef$`Pr(>|t|)`[3]
    p.cis=coef$`Pr(>|t|)`[2]
    
    sig.snp[i]=ifelse(p.snp<alpha, 1, 0)
    sig.cis[i]=ifelse(p.cis<alpha, 1, 0)
    
    simu.coefs[,i]=coef$Estimate
    
  }
  
  
  if(plot.it==TRUE){
    
    
    library('ggplot2', lib="/mnt/ceph/jarredk/Rpackages")
    library('ggridges', lib="/mnt/ceph/jarredk/Rpackages")
    
    if(which.mod=="GMAC"){
      rm.c=match(c("(Intercept)", "pcr1", "platform1","sex2"),row.names(simu.coefs))
    }else{
      rm.c=match(c("(Intercept)"),row.names(simu.coefs))
    }

    long.data=longform(iters = iters, 
                       varvec1 = row.names(simu.coefs)[-c(rm.c)], 
                       cutdata = simu.coefs[-c(rm.c),])
    
    act.data=cbind.data.frame(Coef=names(b)[-rm.c], Value=b[-rm.c])
    
    if(which.mod=="GMAC"){
      fn1=paste0("/mnt/ceph/jarredk/GMACanalysis/permreg_plots/GMAC.ridgeline.trio.", save.name, ".png")
      fn2=paste0("/mnt/ceph/jarredk/GMACanalysis/permreg_plots/GMAC.boxplot.trio.", save.name, ".png")
    }else{
      fn1=paste0("/mnt/ceph/jarredk/GMACanalysis/permreg_plots/ADDIS.ridgeline.trio.", save.name, ".png")
      fn2=paste0("/mnt/ceph/jarredk/GMACanalysis/permreg_plots/ADDIS.boxplot.trio.", save.name, ".png")
    }
    
    png(fn1)
    
    p=ggplot(long.data, aes(x = Value, y = Coef)) +
      geom_density_ridges() +
      geom_point(data=act.data, mapping=aes(x = Value, y = Coef, color = "red"))+
      theme_ridges() + 
      #scale_x_continuous(breaks = seq(-0.2,0.2,0.01), lim = c(-0.2, 0.2))+
      theme(legend.position = "none")
    
    print(p)
    
    dev.off()
    
    
    
    png(fn2)
    
    p=ggplot(long.data, aes(x = Value, y = Coef,)) + 
      geom_boxplot() +
      geom_point(data=act.data, mapping=aes(x = Value, y = Coef, color = "red"))+
      #scale_x_continuous(breaks = seq(-0.2,0.2,0.01), lim = c(-0.2, 0.2))+
      
      theme(legend.position = "none")
    
    print(p)
    
    dev.off()
    
    
    
  }
  
  sig.ind=cbind.data.frame(sig.snp, sig.cis)
  
  return(list(sig.indicator.mat=sig.ind, simu.coefs=simu.coefs, act.coefs=b))
  
}
























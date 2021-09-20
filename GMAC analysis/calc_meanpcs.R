
#calculation of the mean number of adaptively selected pcs for each tissue for all trios
sd.num.pc=matrix(0, nrow = 5, ncol = 4)
colnames(sd.num.pc)=c("unm.cis", "m.cis","unm.trans","m.trans")
row.names(sd.num.pc)=tissues.vec[,1]
mean.num.pc=matrix(0, nrow=5, ncol = 4)
colnames(mean.num.pc)=c("unm.cis", "m.cis","unm.trans","m.trans")
row.names(mean.num.pc)=tissues.vec[,1]

for(i in 1:length(tissues.vec[,1])){
  
  #read in GMAC results
  output.cis=loadRData(fileName=paste('/mnt/ceph/jarredk/GMACanalysis/', tissues.vec[i,1], '/all_trios_output_cis.Rdata', sep = ""))
  output.trans=loadRData(fileName=paste('/mnt/ceph/jarredk/GMACanalysis/', tissues.vec[i,1], '/all_trios_output_trans.Rdata', sep=""))
  
  
  l1=cross.analyze(tissues=tissues.vec[i,1], which.type="cis", save=TRUE)

  l2=cross.analyze(tissues=tissues.vec[i,1], which.type="trans", save=TRUE)
  
  m.cis=rowSums(output.cis$cov.indicator.list[l1$fl.m[[1]]$index.in.trio.mat, ])
  unm.cis=rowSums(output.cis$cov.indicator.list[l1$fl.unm[[1]]$index.in.trio.mat, ])
  m.trans=rowSums(output.trans$cov.indicator.list[l2$fl.m[[1]]$index.in.trio.mat, ])
  unm.trans=rowSums(output.trans$cov.indicator.list[l2$fl.unm[[1]]$index.in.trio.mat, ])
  
  sd.num.pc[i,1]=sd(unm.cis)
  sd.num.pc[i,2]=sd(m.cis)
  sd.num.pc[i,3]=sd(unm.trans)
  sd.num.pc[i,4]=sd(m.trans)
  
  mean.num.pc[i,1]=mean(unm.cis)
  mean.num.pc[i,2]=mean(m.cis)
  mean.num.pc[i,3]=mean(unm.trans)
  mean.num.pc[i,4]=mean(m.trans)
  
  
}
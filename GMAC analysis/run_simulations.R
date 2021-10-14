





#---------------------simulation-tables:-------------------------

source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")

data.tables=vector("list", length = length(tissues.vec[,1]))
nn=5

#small truth model large inferred model
for(i in 1:length(tissues.vec[,1])){
  
  l1=cross.analyze(tissues=tissues.vec[i,1], save=FALSE)
  all.m0m3=l1$final.tables[[1]][c(which(l1$final.tables[[1]]$Addis.Class=="M3"), which(l1$final.tables[[1]]$Addis.Class=="MO")),]
  print(head(all.m0m3[,1:6]))
  trios=sample(all.m0m3$Trio.Num, nn)
  ot1=run.simu12(tissue = tissues.vec[i,1] ,trios=trios, 
                 mod.type.vec=all.m0m3$Mediation.type[match(trios, all.m0m3$Trio.Num)],
                 alpha=0.001, n=10)
  ot1$out.mat$Tissue=rep(tissues.vec[i,1])
  data.tables[[i]]=ot1$out.mat
  
}


final.table=apply(data.tables, 1, rbind.data.frame)
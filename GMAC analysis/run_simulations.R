





#---------------------simulation-tables:-------------------------

source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")


nn=50

#small truth model large inferred model
for(i in 1:length(tissues.vec[,1])){
  
  print(tissues.vec[i,1])
  l1=cross.analyze(tissues=tissues.vec[i,1], save=FALSE)
  all.m0m3=l1$final.tables[[1]][c(which(l1$final.tables[[1]]$Addis.Class=="M3"), which(l1$final.tables[[1]]$Addis.Class=="MO")),]
  print(head(all.m0m3[,1:6]))
  trios=sample(all.m0m3$Trio.Num, nn)
  print("running Simulations...")
  ot1=run.simu12(tissue = tissues.vec[i,1] ,trios=trios, 
                 l1.table=all.m0m3,
                 mod.type.vec=all.m0m3$Mediation.type[match(trios, all.m0m3$Trio.Num)],
                 alpha=0.001, n="random")
  print("...done")
  ot1$Tissue=rep(tissues.vec[i,1])
  
  if(i == 1){
    write.table(ot1, file = "/mnt/ceph/jarredk/GMACanalysis/master_tables/simulations.txt", col.names = T, row.names = F, sep = "\t", 
                quote=F)
  }else{
    write.table(ot1, file = "/mnt/ceph/jarredk/GMACanalysis/master_tables/simulations.txt", col.names = F, row.names = F, sep = "\t",
                quote=F, append = T)
  }
  
}



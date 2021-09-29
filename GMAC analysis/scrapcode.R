
  #read in postprocessing functions

source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
l1=cross.analyze(tissues=tissues.vec[,1], save=FALSE)
l1$final.tables[[5]][which(l1$final.tables[[5]]$Addis.Class=="M3"),][1:10,]
Lond2Addis.lookup(trio.index=135, tissue.name="WholeBlood", with.pc=TRUE)[2:4]
ld1=runit(trio.number = 133, mtype = "M3", pcs=paste0("PC",c(11,21,24)) )
ld1=runit(trio.number = 135, mtype = "M3", pcs=NULL )

pcor(ld1$GMAC)$estimate[1:2,1:3]
pcor(ld1$addis)$estimate[1:2,1:3]


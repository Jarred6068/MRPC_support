#read in postprocessing functions
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
l1=cross.analyze(tissues=tissues.vec[5,1], save=FALSE)

#check first 15 M0's for wholeblood
data2=l1$final.tables[[1]][which(l1$final.tables[[1]]$Addis.Class=="MO"),]
M0.trios.sample=l1$final.tables[[1]][which(l1$final.tables[[1]]$Addis.Class=="MO"),]$Trio.Num

x=c(2697,2914,3541,4747,5195,5822,5919,6669,6756)

prb.ind=match(x, M0.trios.sample)

outM0=runit(indata = data2[-prb.ind,], trio.number =  M0.trios.sample[-prb.ind], mtype = "M0")

write.csv(outM0$table, file = "/mnt/ceph/jarredk/GMACanalysis/master_tables/trios.pcors.summary.tableM0.csv", quote = FALSE,
          row.names = FALSE, col.names = TRUE)

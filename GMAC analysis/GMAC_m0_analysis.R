#read in postprocessing functions
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
l1=cross.analyze(tissues=tissues.vec[5,1], save=FALSE)

#check first 15 M0's for wholeblood
data2=l1$final.tables[[1]][which(l1$final.tables[[1]]$Addis.Class=="MO"),]
M0.trios.sample=l1$final.tables[[1]][which(l1$final.tables[[1]]$Addis.Class=="MO"),]$Trio.Num

outM0=runit(indata = data2[1:30,], trio.number =  M0.trios.sample[1:30], mtype = "M0")

write.csv(outM0$table, file = "/mnt/ceph/jarredk/GMACanalysis/master_tables/trios.pcors.summary.tableM0.csv", quote = FALSE,
          row.names = FALSE, col.names = TRUE)

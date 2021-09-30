#read in postprocessing functions
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
l1=cross.analyze(tissues=tissues.vec[5,1], save=FALSE)
#check first 15 M3's for wholeblood
data1=l1$final.tables[[1]][which(l1$final.tables[[1]]$Addis.Class=="M3"),]
M3.trios.sample=l1$final.tables[[1]][which(l1$final.tables[[1]]$Addis.Class=="M3"),]$Trio.Num


outM3=runit(indata = data1[1:30,], trio.number =  M3.trios.sample[1:30], mtype = "M3")

write.csv(outM3$table, file = "/mnt/ceph/jarredk/GMACanalysis/master_tables/trios.pcors.summary.tableM3.csv", quote = FALSE,
          row.names = FALSE, col.names = TRUE)


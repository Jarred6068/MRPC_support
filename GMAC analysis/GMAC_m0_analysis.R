#read in postprocessing functions
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
l1=cross.analyze(tissues=tissues.vec[5,1], save=FALSE)

#check first 15 M0's for wholeblood
data2=l1$final.tables[[1]][which(l1$final.tables[[1]]$Addis.Class=="MO"),]


outM0=runit(indata = data2, trio.number =  data2$Trio.Num, med.type=data2$Mediation.type, mtype = "M0")

write.csv(outM0$table, file = "/mnt/ceph/jarredk/GMACanalysis/master_tables/trios.pcors.summary.tableM0.csv", quote = FALSE,
          row.names = FALSE, col.names = TRUE)

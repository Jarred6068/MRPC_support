#read in postprocessing functions
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
l1=cross.analyze(tissues=tissues.vec[5,1], save=FALSE)
#check first 15 M3's for wholeblood
data1=l1$final.tables[[1]][which(l1$final.tables[[1]]$Addis.Class=="M3"),]
M3.trios.sample=l1$final.tables[[1]][which(l1$final.tables[[1]]$Addis.Class=="M3"),]$Trio.Num

#check first 15 M0's for wholeblood
data2=l1$final.tables[[1]][which(l1$final.tables[[1]]$Addis.Class=="M0"),]
M0.trios.sample=l1$final.tables[[1]][which(l1$final.tables[[1]]$Addis.Class=="M0"),]$Trio.Num

#Lond2Addis.lookup(trio.index=135, tissue.name="WholeBlood", with.pc=TRUE)[2:4]
outM3=runit(indata = data1[c(1,3,5:12),], trio.number =  M3.trios.sample[c(1,3,5:12)], mtype = "M3")
#outM0=runit(indata = data2[1:15,], trio.number =  M0.trios.sample[1:15], mtype = "M0")
#ld1=runit(trio.number = 135, mtype = "M3", pcs=NULL )

write.csv(outM3$table, file = "/mnt/ceph/jarredk/GMACanalysis/trios.pcors.summary.tableM3.csv", quote = FALSE,
          row.names = FALSE, col.names = TRUE)

#write.csv(outM0$table, file = "/mnt/ceph/jarredk/GMACanalysis/trios.pcors.summary.tableM0.csv", quote = FALSE,
#          row.names = FALSE, col.names = TRUE)

#pcor(ld1$GMAC)$estimate[1:2,1:3]
#pcor(ld1$addis)$estimate[1:2,1:3]

#write master table to csv
# final.table=rbind.data.frame(outM3$table, outM0$table)
# 
# write.csv(final.table, file = "/mnt/ceph/jarredk/GMACanalysis/30trios.pcors.summary.table.csv", quote = FALSE,
#           row.names = FALSE, col.names = TRUE)



source("/mnt/ceph/jarredk/Reg_Net/AL_genetabV2.R")

sort.gt2(tissue.selection=c(1:48), FDR="ADDIS", mediator="cis", save.data=TRUE, verbose=TRUE)
sort.gt2(tissue.selection=c(1:48), FDR="ADDIS", mediator="trans", save.data=TRUE, verbose=TRUE)
sort.gt2(tissue.selection=c(1:48), FDR="LOND", mediator="cis", save.data=TRUE, verbose=TRUE)
sort.gt2(tissue.selection=c(1:48), FDR="LOND", mediator="trans", save.data=TRUE, verbose=TRUE)
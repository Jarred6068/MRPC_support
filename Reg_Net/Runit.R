

source("/mnt/ceph/jarredk/Reg_Net/AL_genetabV2.R")
sort.gt2()
sort.gt2(mediator="cis")
sort.gt2(FDR="LOND",mediator="cis")
sort.gt2(FDR="LOND",mediator="trans")

df=master.create()
dim(df$lm1t1)
dim(df$lm1t2)
dim(df$am1t1)
dim(df$am1t2)

trans.as.mediator=binit(df$am1t2, target=4)
cis.as.mediator=binit(df$am2t1, target=3)



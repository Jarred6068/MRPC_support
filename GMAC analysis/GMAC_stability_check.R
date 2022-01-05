
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")

output.WB.sc=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/WholeBlood/all_trios_output_cisstability_check.Rdata')
output.AS.sc=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/AdiposeSubcutaneous/all_trios_output_cisstability_check.Rdata')
output.t.WB.sc=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/WholeBlood/all_trios_output_transstability_check.Rdata')
output.t.AS.sc=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/AdiposeSubcutaneous/all_trios_output_trans_stability_check.Rdata')

WB.sc=output.WB.sc$output.table$pval_Known_sel_pool
AS.sc=output.AS.sc$output.table$pval_Known_sel_pool

WB.t.sc=output.t.WB.sc$output.table$pval_Known_sel_pool
AS.t.sc=output.t.AS.sc$output.table$pval_Known_sel_pool

l1=cross.analyze(tissues=c("WholeBlood","AdiposeSubcutaneous"), save=FALSE)$final.tables

#pull out the first run pvalues 
#cis trios
WB=l1[[1]][which(l1[[1]]$Mediation.type!="Trans.Med"),]
AS=l1[[2]][which(l1[[1]]$Mediation.type!="Trans.Med"),]
#trans trios
WB.t=l1[[1]][which(l1[[1]]$Mediation.type=="Trans.Med"),]
AS.t=l1[[2]][which(l1[[1]]$Mediation.type=="Trans.Med"),]

#check cis trios
row.match(output.WB.sc$output.table[,1:3], WB[,1:3])
comb.WB=cbind.data.frame(WB.sc, WB$Perm.reg.p, WB$Trio.Num)
comb.WB[which(comb.WB[,1]>0.05),]
#check trans trios
row.match(output.t.WB.sc$output.table[,1:3], WB.t[,1:3])
comb.WB.t=cbind.data.frame(WB.t.sc, WB.t$Perm.reg.p, WB.t$Trio.Num)
comb.WB.t[which(comb.WB.t[,1]>0.05),]


row.match(output.AS.sc$output.table[,1:3], AS[,1:3])
comb.AS=cbind.data.frame(AS.sc, AS$Perm.reg.p, AS$Trio.Num)
comb.AS[which(comb.AS[,1]>0.05),]


trios1=c(2417, 8662, 4268, 8112)
trios2=c(1930, 2573, 3374, 7770, 8308)





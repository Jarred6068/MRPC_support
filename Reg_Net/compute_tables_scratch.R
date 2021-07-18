
#import the deprecated genes:

miss1=read.csv(file="/mnt/ceph/jarredk/Reg_Net/list_of_deprecated_genes.csv")

#combine ensembl 37 and 38 biomart datasets by the adding to 38 the genes unique to 37
meta.data=read.table(file="/mnt/ceph/jarredk/Reg_Net/mart_export37.txt", sep="\t", header=T)
meta.data2=read.table(file="/mnt/ceph/jarredk/Reg_Net/mart_export38.txt", sep="\t", header=T)
nomatch=which(is.na(match(meta.data$Gene.stable.ID, meta.data2$Gene.stable.ID)))
merged.meta.data=rbind.data.frame( miss1[,-1], meta.data[nomatch,], meta.data2)

write.table(merged.meta.data, file = "/mnt/ceph/jarredk/Reg_Net/mart_export_merged.txt", sep = "\t", quote = F, row.names = F, col.names = T)



#write all tissue tables
source("/mnt/ceph/jarredk/Reg_Net/AL_genetabV2.R")
sort.gt2(tissue.selection=c(1:48), FDR="ADDIS", mediator="trans", save.data=FALSE, verbose=TRUE)
sort.gt2(tissue.selection=c(1:48), FDR="ADDIS", mediator="cis", save.data=FALSE, verbose=TRUE)
sort.gt2(tissue.selection=c(1:48), FDR="LOND", mediator="trans", save.data=FALSE, verbose=TRUE)
sort.gt2(tissue.selection=c(1:48), FDR="LOND", mediator="cis", save.data=FALSE, verbose=TRUE)

df=master.create()

save(df, file = "/mnt/ceph/jarredk/Reg_Net/loaded_AL_datatables.Rdata")








#read in healthy breast tissue methylation data
hbt.data=read.table(file="/mnt/ceph/jarredk/Cancer_Selection/GSE101961_normal_breast/GSE101961.unnorm.methyl.final.txt", header=T, sep="\t")
hbt.pvals=read.table(file="/mnt/ceph/jarredk/Cancer_Selection/GSE101961_normal_breast/GSE101961.signal.pvals.txt", header=T, sep="\t")
hbt.meta=read.csv(file="/mnt/ceph/jarredk/Cancer_Selection/GSE101961_normal_breast/GSE101961_meta_info.csv")
dim(hbt.data)
#read in stem_cell data set 1
stem.meth1=read.table("/mnt/ceph/jarredk/Cancer_Selection/GSE116754_METH_stem/GSE59091/GSE59091_prop_methyl.txt", sep="\t", header=T)
stem.meth1.meta=read.csv("/mnt/ceph/jarredk/Cancer_Selection/GSE116754_METH_stem/GSE59091/GSE59091_meta_info.csv")
dim(stem.meth1)
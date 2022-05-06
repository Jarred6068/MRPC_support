
#specify package path:
pack.path="/mnt/ceph/jarredk/Rpackages"

#read in all normalized data sets:
#read in healthy breast tissue methylation data
hbt.data=read.table("/mnt/ceph/jarredk/Cancer_Selection/GSE101961_normal_breast/GSE101961_normal_breast_log_normal.txt",
                    header = T, sep = "\t")
#hbt.data=read.table("/mnt/ceph/jarredk/Cancer_Selection/GSE101961_normal_breast/GSE101961_normal_breast_log_normal.txt",
#                    header = T, sep = "\t")
hbt.labels=rep("Normal", dim(hbt.data)[2])

#read in stem_cell data set 59091
stc59=read.table("/mnt/ceph/jarredk/Cancer_Selection/GSE116754_METH_stem/GSE59091/GSE59091_log_normal.txt",
                 header = T, sep = "\t")
#stc59=read.table("/mnt/ceph/jarredk/Cancer_Selection/GSE116754_METH_stem/GSE59091/GSE59091_log_normal.txt",
#                 header = T, sep = "\t")
stc59.meta=read.csv("/mnt/ceph/jarredk/Cancer_Selection/GSE116754_METH_stem/GSE59091//GSE59091_meta_info.csv")

#read in stem_cell data set 116754
#stc11=read.table("/mnt/ceph/jarredk/Cancer_Selection/GSE116754_METH_stem/GSE116754/GSE116754_STEM_log_normal.txt",
#                 header = T, sep = "\t")
stc11=read.table("/mnt/ceph/jarredk/Cancer_Selection/GSE116754_METH_stem/GSE116754/GSE116754_STEM_log_normal.txt",
                 header = T, sep = "\t")
stc11.meta=read.csv("/mnt/ceph/jarredk/Cancer_Selection/GSE116754_METH_stem/GSE116754/GSE116754_meta_info.csv")

#breast cancer data
bc.data=read.table("/mnt/ceph/jarredk/Cancer_Selection/GSE75067_BC_with_annotations/GSE75067_BC_log_normal.txt",
                   header = T, sep = "\t")
#bc.data=read.table("/mnt/ceph/jarredk/Cancer_Selection/GSE75067_BC_with_annotations/GSE75067_BC_log_normal.txt",
#                   header = T, sep = "\t")
bc.meta=read.table("/mnt/ceph/jarredk/Cancer_Selection/GSE75067_BC_with_annotations/GSE75067_sample_annotations.txt",
                   header = T, sep = "\t")
#get the list of samples to remove for having missing grade or ER status
bc.grade.ER=cbind.data.frame(ER.status=bc.meta$ER, hist.grade=bc.meta$grade)
samples.rm.bc=unique(attr(na.omit(bc.grade.ER),"na.action"))
bc.grade.ER[samples.rm.bc,]
table(bc.grade.ER$ER.status, bc.grade.ER$hist.grade)

bc.labels=paste0(bc.meta$ER, bc.meta$grad)[-samples.rm.bc]

#align all data.sets and combine into large df:
#smallest is stc59
#############################################
align.set=function(input=NULL, target=NULL){
  idx=match(row.names(input), row.names(target))
  X=input[na.omit(idx),]
  return(X)
}
#############################################

bc.data.aligned=align.set(input = bc.data, target = stc59)
hbt.data.aligned=align.set(input = hbt.data, target = stc59)
stc11.aligned=align.set(input = stc11, target = stc59)
#keep only HESc samples in stc11 and non-reprogrammed cells in stc59:
samples.kept.stc59=which(stc59.meta$Reprogram_Method==" NA")
samples.kept.stc11=which(stc11.meta$Tissue==" Human Embryonic Stem Cells")
#combine
all.data=cbind.data.frame(stc59[,samples.kept.stc59],
                          stc11.aligned[,samples.kept.stc11],
                          hbt.data.aligned,
                          bc.data.aligned[,-samples.rm.bc])

#save(all.data, file = "/mnt/ceph/jarredk/Cancer_Selection/all.data.final.RData")

labels.final=c(stc59.meta$Donor_Cell_Type[samples.kept.stc59],
               stc11.meta$Tissue[samples.kept.stc11],
               hbt.labels,
               bc.labels)
#fix labels
labels.final[which(labels.final==" Human Embryonic Stem Cells")]="HESC"
labels.final[which(labels.final==" embryonic stem cell")]="HESC"
labels.final[which(labels.final==" dermal fibroblast")]="dermal.fibroblast"
labels.final[which(labels.final==" endothelial precursor")]="endothelial.precursor"
labels.final[which(labels.final==" foreskin fibroblast")]="foreskin.fibroblast"

colnames(all.data)=labels.final
save(all.data, file = "/mnt/ceph/jarredk/Cancer_Selection/in_depth_analysis/all.data.final.prop.methyl.RData")

save(labels.final, file = "/mnt/ceph/jarredk/Cancer_Selection/in_depth_analysis/all.data.final.labels.prop.methyl.RData")

#
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager", lib=pack.path)
#
# BiocManager::install("IlluminaHumanMethylation450k.db", lib=pack.path)
#
# library("IlluminaHumanMethylation450k.db", lib=pack.path)

#specify package path:
pack.path="/mnt/ceph/jarredk/Rpackages"
#load libraries
library('NbClust', lib = pack.path)
#library('readxl', lib = pack.path)
library('qvalue', lib = pack.path)
#read in human methylation 450k header descriptions
hm450=read.csv("/mnt/ceph/jarredk/Cancer_Selection/GPL13534_HumanMethylation450_15017482_v.1.1.csv")

load(file = "/mnt/ceph/jarredk/Cancer_Selection/in_depth_analysis/all.data.final.prop.methyl.RData")

load(file = "/mnt/ceph/jarredk/Cancer_Selection/in_depth_analysis/all.data.final.labels.prop.methyl.RData")

genes=hm450$UCSC_RefGene_Name[match(row.names(all.data), hm450$IlmnID)]
genes.list=strsplit(genes, ";")
#replace character(0) with NA
genes.list.NAs=lapply(genes.list, function(x){ifelse(identical(x, character(0)), NA, x)})
nogenenames=which(is.na(genes.list.NAs))

all.data=all.data[-nogenenames,]
genes.shortlist=genes.list[-nogenenames]
unq.genes=as.list(unique(unlist(genes.shortlist)))

#############################################
get.gene.idx=function(unique.genes=NULL, genes.list=NULL){
  gene.idx.list=list()
  for(i in 1:length(unq.genes)){
    gene.idx.list[[i]]=which(lapply(genes.list, function(x,y){any(match(x, y))}, y=unique.genes[i])==TRUE)
    #print(gene.idx.list[[i]])
  }

  return(gene.idx.list)

}
#############################################

mean.gene=function(idx=NULL, full.data=NULL){

  if(length(idx)>1){
    y=colMeans(full.data[idx,], na.rm = TRUE)
  }else{
    y=full.data[idx,]
  }


}

#############################################

calc.centroids=function(gene.idx.list=NULL, unique.genes=NULL, full.data=NULL){

  gene.centroids=sapply(gene.idx.list, mean.gene, full.data=full.data, simplify = "array")
  print(dim(gene.centroids))
  V=matrix(unlist(gene.centroids), nrow=333, ncol=21244, byrow=F)
  colnames(V)=unique.genes
  row.names(V)=colnames(full.data)
  return(V)
}

#############################################
#get gene.idx list
#unq.gene.rowid=get.gene.idx(unique.genes = unq.genes, genes.list = genes.shortlist)
#save(unq.gene.rowid, file="/mnt/ceph/jarredk/Cancer_Selection/unq.gene.rowid.RData")
load(file="/mnt/ceph/jarredk/Cancer_Selection/unq.gene.rowid.RData")
#get the centroids matrix
all.data.centroids=calc.centroids(gene.idx.list = unq.gene.rowid, unique.genes = unq.genes, full.data = all.data)
save(all.data.centroids, file = "/mnt/ceph/jarredk/Cancer_Selection/in_depth_analysis/all.data.centroids2.RData")

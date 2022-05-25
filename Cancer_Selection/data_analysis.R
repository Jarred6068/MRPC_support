#specify package path:
pack.path="/mnt/ceph/jarredk/Rpackages"
#load libraries
library('NbClust', lib = pack.path)
#library('readxl', lib = pack.path)
library('qvalue', lib = pack.path)
#read in human methylation 450k header descriptions
hm450=read.csv("/mnt/ceph/jarredk/Cancer_Selection/GPL13534_HumanMethylation450_15017482_v.1.1.csv")

load(file = "/mnt/ceph/jarredk/Cancer_Selection/all.data.final.RData")

load(file = "/mnt/ceph/jarredk/Cancer_Selection/all.data.final.labels.RData")

#remove probe location variable
sites.on.genome=strsplit(hm450$UCSC_RefGene_Group,";")
#get TSS sites and gene body site probe indexes
TSS.sites=which(unlist(lapply(sites.on.genome, function(x){any(grepl("TSS", x))})))
gene.body.sites=which(unlist(lapply(sites.on.genome, function(x){any(grepl("Body", x))})))
all.data.new=all.data[,]

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
save(all.data.centroids, file = "/mnt/ceph/jarredk/Cancer_Selection/all.data.centroids.RData")
load(file = "/mnt/ceph/jarredk/Cancer_Selection/all.data.centroids.RData")
all.data.centroids=t(all.data.centroids)

#get differentially methylated genes:
#first remove the stem cell samples:
all.data.diff=all.data.centroids[,-c(1:56)] #first 57 columns are from the stem cell data
labels.diff=labels.final[-c(1:56)]
#labels.diff[labels.diff!="Normal"]="Cancer"

#############################################
get.anova.p=function(X, factor1){

  X.new=cbind.data.frame(gene=X, type=as.factor(factor1))
  #print(head(X.new))
  anova.gene=anova(lm(gene~., data=X.new))
  #print(anova.gene)
  return(anova.gene$`Pr(>F)`[1])

}

#############################################

get.diff.genes=function(centroid.data=NULL, correction=c("qvalue","BH","bonferroni","none"), factor.labels=NULL,
                        FDR=0.05, alpha=0.05){

  X=t(centroid.data)

  p.values=apply(X, 2, get.anova.p, factor1=factor.labels)
  names(p.values)=colnames(X)

  adj.p=NULL
  q.val=NULL

  switch(correction, qvalue = {

    q.val=qvalue(p.values, fdr.level = FDR)$qvalue
    sig.genes=names(q.val[q.val<alpha])

  }, BH = {

    adj.p=p.adjust(p.values, method = correction)
    sig.genes=names(adj.p[adj.p<alpha])

  }, bonferroni = {

    adj.p=p.adjust(p.values, method = correction)
    sig.genes=names(adj.p[adj.p<alpha])

  }, none = {

    sig.genes=names(p.values[p.values<alpha])

  },stop("Correction not included or missing"))

  return(list(raw.p=p.values,
              adjusted.p=adj.p,
              q.values=q.val,
              genes=sig.genes,
              G=X[,which(q.val<alpha)]))

}
#############################################

#run anova and return obj. with pvalues and sig.gene.matrix
out1=get.diff.genes(centroid.data = all.data.diff,
                    correction = 'qvalue',
                    factor.labels = labels.diff,
                    FDR=0.05,
                    alpha=0.00005)

save(out1, file="/mnt/ceph/jarredk/Cancer_Selection/output1.RData")

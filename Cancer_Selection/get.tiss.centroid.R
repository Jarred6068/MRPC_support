
setwd("/mnt/ceph/jarredk/Cancer_Selection")
dm=read.table("differential.matrix.txt",header=T)
load("all.data.final.labels.RData")
tiss=labels.final[-c(1:56)]
x=unique(tiss)
tiss.centroids=as.data.frame(matrix(0, nrow = length(x), ncol = ncol(dm)))
for(i in 1:length(x)){

  row.idx=which(tiss==x[i])
  tiss.centroids[i, ]=colMeans(dm[row.idx,], na.rm = T)
  row.names(tiss.centroids)=x

}

library('ape', lib="/mnt/ceph/jarredk/Rpackages")

tree1=nj(dist(tiss.centroids, method = "euclidean"))
pdf()
plot.phylo(tree1)
dev.off()

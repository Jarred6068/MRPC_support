
setwd("/mnt/ceph/jarredk/Cancer_Selection/Treeplots")
library('ape', lib="/mnt/ceph/jarredk/Rpackages")
library('phytools',lib="/mnt/ceph/jarredk/Rpackages")
dm=read.table("/mnt/ceph/jarredk/Cancer_Selection/differential.matrix.txt",header=T)
load("/mnt/ceph/jarredk/Cancer_Selection/all.data.final.labels.RData")
tiss=labels.final
x=unique(tiss)
tiss.centroids=as.data.frame(matrix(0, nrow = length(x), ncol = ncol(dm)))
for(i in 1:length(x)){

  row.idx=which(tiss==x[i])
  tiss.centroids[i, ]=colMeans(dm[row.idx,], na.rm = T)
  row.names(tiss.centroids)=x

}

save(tiss.centroids, file = "tissue.centroids.diff.genes.RData")



#easy function to wrap several distance metric calculations
calc.dist=function(x = NULL, use.method = c("corr","euclidean","minkowski","maximum")){

  switch(use.method, corr = {

    distmat = 1-cor(t(x))^2
    return(distmat)

  }, euclidean = {

    distmat = dist(x, method = "euclidean")
    return(distmat)

  }, minkowski = {

    distmat = dist(x, method = "minkowski", p=1)
    return(distmat)

  }, maximum = {

    distmat = dist(x, method = "maximum")
    return(distmat)

  })


}


#function to calculate WLS and NJ trees for each distance
build.trees=function(x=NULL){
  trees=list()
  dist.types=c("corr","euclidean","minkowski","maximum")
  for(i in 1:4){

    nj.trees=list()
    me.trees=list()
    for(j in 1:10000){
      sampmat = x[,sample(c(1:ncol(x)), round(0.7*ncol(x)))]
      D=calc.dist(sampmat, use.method = dist.types[i])
      nj.trees[[j]]=nj(as.dist(D))
      me.trees[[j]]=fastme.ols(as.dist(D))
    }

    save(nj.trees, file = paste0("NJ.",dist.types[i], ".RData"))
    save(me.trees, file = paste0("ME_WLS..",dist.types[i], ".RData"))
    #nj.ctree=compute.brlen(consensus(nj.trees, rooted = F, p=0.95))
    #me.ctree=compute.brlen(consensus(me.trees, rooted = F, p=0.95))
    nj.ctree=consensus.edges(nj.trees, method="mean.edge")
    me.ctree=consensus.edges(me.trees, method="mean.edge")
    print(nj.ctree)
    str(nj.ctree)
    root.nj.ctree=root(nj.ctree, outgroup="HESC", resolve.root=T)
    root.me.ctree=root(me.ctree, outgroup="HESC", resolve.root=T)
    print(root.nj.ctree$edge.length)
    print(root.me.ctree$edge.length)
    print(is.rooted(root.nj.ctree))
    print(is.rooted(root.me.ctree))

    print(paste0("Plotting concensus trees using metric ", dist.types[i]))
    pdf(paste0("nj.me.ctree.nostem.",dist.types[i],".pdf"))
    par(mfrow=c(1,2))
    plot.phylo(root.nj.ctree,
               main = paste0("Neighbor-Joining ", dist.types[i]), root.edge = T)
    plot.phylo(root.me.ctree,
               main = paste0("Minimum Evolution WLS ",dist.types[i]), root.edge = T)
    par(mfrow=c(1,1))
    dev.off()

    trees[[i]]=list(nj=root.nj.ctree, me=root.me.ctree)
  }

  return(trees)
}

tl=build.trees(x=tiss.centroids[-c(2:4),])


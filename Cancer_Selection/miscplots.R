load("/mnt/ceph/jarredk/Cancer_Selection/output1.RData")
load("/mnt/ceph/jarredk/Cancer_Selection/all.data.centroids.RData")
load("/mnt/ceph/jarredk/Cancer_Selection/all.data.final.RData")
load("/mnt/ceph/jarredk/Cancer_Selection/all.data.final.labels.RData")

library('ggpubr', lib="/mnt/ceph/jarredk/Rpackages")


pdf("/mnt/ceph/jarredk/Cancer_Selection/raw.p.hist.pdf")
gghistogram(out1$raw.p)
dev.off()




data.diff=all.data.centroids[-c(1:56),]
diff.labels=labels.final[-c(1:56)]
ss=sample(c(1:dim(data.diff)[1]), 100)
for(i in 1:100){
  pdf(paste0("/mnt/ceph/jarredk/Cancer_Selection/Boxplots/bp",i,".pdf"))
  data1=cbind.data.frame(labels = diff.labels, gene=data.diff[,ss[i]])
  boxplot(gene~labels, data=data1, main = colnames(data.diff)[ss[i]],
          xlab = "Tissue Condition",
          ylab = "Methylation")
  dev.off()
}



#run anova and return obj. with pvalues and sig.gene.matrix
out1=get.diff.genes(centroid.data = t(data.diff),
                    correction = 'qvalue',
                    factor.labels = diff.labels,
                    FDR=0.05,
                    alpha=0.01)

library(ggpubr, lib="/mnt/ceph/jarredk/Rpackages")
bmart=read.table(file="/mnt/ceph/jarredk/Reg_Net/mart_export_merged_lncRNA_fixed.txt", sep = "\t", header = T)
setwd("/mnt/ceph/jarredk/Cancer_Selection")

dm=read.table("/mnt/ceph/jarredk/Cancer_Selection/differential.matrix.txt",header=T)
load("/mnt/ceph/jarredk/Cancer_Selection/all.data.final.labels.RData")
row.names(dm)=make.unique(labels.final)
tissue=as.factor(labels.final)
pca=prcomp(dm[-c(1:56),], retx=TRUE)
rot.var=as.data.frame(pca$x)
rn=row.names(pca$rotation)
pc1=pca$rotation[,1]
pos.genes=rn[which(pc1<0)]
neg.genes=rn[which(pc1>0)]
types.pos=bmart$Gene.type[match(pos.genes, bmart$Gene.name)]
types.neg=bmart$Gene.type[match(neg.genes, bmart$Gene.name)]
#print(pos.genes)
#print(neg.genes)
ev=as.data.frame(pca$rotation)
row.names(ev)=rn
break.down.pos=summary(as.factor(types.pos))
break.down.neg=summary(as.factor(types.neg))

types.pc1.pl=cbind.data.frame(Loading.on.PC1=c(rep("pos.loading", length(break.down.pos)),
                                               rep("neg.loading", length(break.down.neg))),
                              gene.types=c(names(break.down.pos), names(break.down.neg)),
                              counts=c(log10(break.down.pos), log10(break.down.neg)))

v=pca$sdev^2/sum(pca$sdev^2)
dat2=cbind.data.frame(Percent.of.Variance=v[1:9],
                      Principle.Component=paste0("PC",c(1:9)))



#show the log.p scale
p=sort(runif(1000))
log.p=log(p/(1-p))
G=cbind.data.frame(p=p, log.p=log.p)
pdf("/mnt/ceph/jarredk/Cancer_Selection/in_depth_analysis/logpkey.pdf")
A=ggplot(data=G,aes(x=p, y=log.p))+
  geom_line(size=3)+
  xlab("Proportion of Methylation")+
  ylab("Logit(Proportion of Methylation")+
  ggtitle("Logit Transformation")+
  theme_pubr()
plot(A)
dev.off()


rot.var2=cbind.data.frame(rot.var, Hist.Group=tissue[-c(1:56)])
pdf("/mnt/ceph/jarredk/Cancer_Selection/in_depth_analysis/PCA_PC1_2.pdf")
sp <- ggscatter(rot.var2, x = "PC1", y = "PC2",
                color = "Hist.Group", palette = "jco",
                size = 2, alpha = 0.6)+
  border()
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(rot.var2, "PC1", fill = "Hist.Group",
                   palette = "jco")
yplot <- ggdensity(rot.var2, "PC2", fill = "Hist.Group",
                   palette = "jco")+
  rotate()
# Cleaning the plots
yplot <- yplot + clean_theme()
xplot <- xplot + clean_theme()
# Arranging the plot
ggarrange(xplot, NULL, sp, yplot,
          ncol = 2, nrow = 2,  align = "hv",
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)

dev.off()



pdf("/mnt/ceph/jarredk/Cancer_Selection/in_depth_analysis/PCA_PC2_3.pdf")
sp <- ggscatter(rot.var2, x = "PC2", y = "PC3",
                color = "Hist.Group", palette = "jco",
                size = 2, alpha = 0.6)+
  border()
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(rot.var2, "PC2", fill = "Hist.Group",
                   palette = "jco")
yplot <- ggdensity(rot.var2, "PC3", fill = "Hist.Group",
                   palette = "jco")+
  rotate()
# Cleaning the plots
yplot <- yplot + clean_theme()
xplot <- xplot + clean_theme()
# Arranging the plot
ggarrange(xplot, NULL, sp, yplot,
          ncol = 2, nrow = 2,  align = "hv",
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)

dev.off()




#subsidiary information about genes
dat.new=cbind.data.frame(Hist.group=as.factor(labels.final[-c(1:56)]), dm[-c(1:56),])
pdf("/mnt/ceph/jarredk/Cancer_Selection/in_depth_analysis/loading.gene.content.pdf")
C=ggplot(data=types.pc1.pl, aes(x=gene.types, y=counts, fill=Loading.on.PC1)) +
  geom_bar(stat="identity",color="black")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylab("log10(count)")
D=ggboxplot(data=dat.new, x="Hist.group", y="BCL8")+
  scale_fill_grey() + theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


E=ggboxplot(data=dat.new, x="Hist.group", y="GAGE1")+
  scale_fill_grey() + theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggarrange(C,
          ggarrange(D, E, ncol = 2, labels = c("B", "C")),
          nrow = 2,
          labels = "A")
dev.off()

pdf("/mnt/ceph/jarredk/Cancer_Selection/in_depth_analysis/screePCA.pdf")
B=ggplot(data=data.frame(dat2), aes(x=Principle.Component,y=Percent.of.Variance))+
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=round(Percent.of.Variance,3)), vjust=1.2, color="white", size=3.5)+
  theme_minimal()
plot(B)
dev.off()








#LDA genes
Hist.Group2=as.factor(labels.final[-c(1:56)])
library(MASS)
lda.obj=lda(data.frame(dm[-c(1:56),]), grouping=Hist.Group2)

scores=cbind.data.frame(as.data.frame(as.matrix(dm[-c(1:56),])%*%lda.obj$scaling), Histology=Hist.Group2)


pdf("/mnt/ceph/jarredk/Cancer_Selection/in_depth_analysis/LDA_LD1_2.pdf")
sp <- ggscatter(scores, x = "LD1", y = "LD2",
                color = "Histology", palette = "jco",
                size = 2, alpha = 0.6)+
  border()
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(scores, "LD1", fill = "Histology",
                   palette = "jco")
yplot <- ggdensity(scores, "LD2", fill = "Histology",
                   palette = "jco")+
  rotate()
# Cleaning the plots
yplot <- yplot + clean_theme()
xplot <- xplot + clean_theme()
# Arranging the plot
ggarrange(xplot, NULL, sp, yplot,
          ncol = 2, nrow = 2,  align = "hv",
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)

dev.off()























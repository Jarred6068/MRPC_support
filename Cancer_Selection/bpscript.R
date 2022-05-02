#subsidiary information about genes
load("/mnt/ceph/jarredk/Cancer_Selection/in_depth_analysis/all.data.centroids2.RData")
load("/mnt/ceph/jarredk/Cancer_Selection/in_depth_analysis/all.data.final.labels.prop.methyl.RData")
dat.new=cbind.data.frame(Hist.group=as.factor(labels.final), all.data.centroids)
pdf("/mnt/ceph/jarredk/Cancer_Selection/in_depth_analysis/MIR193B_MIR367.pdf")

D=ggboxplot(data=dat.new, x="Hist.group", y="MIR193B")+
  scale_fill_grey() + theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


E=ggboxplot(data=dat.new, x="Hist.group", y="MIR367")+
  scale_fill_grey() + theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggarrange(D, E, ncol = 2, labels = c("A", "B"))
dev.off()



pdf("/mnt/ceph/jarredk/Cancer_Selection/in_depth_analysis/BCL8_GAGE1.pdf")

D=ggboxplot(data=dat.new, x="Hist.group", y="BCL8")+
  scale_fill_grey() + theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


E=ggboxplot(data=dat.new, x="Hist.group", y="GAGE1")+
  scale_fill_grey() + theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggarrange(D, E, ncol = 2, labels = c("A", "B"))
dev.off()

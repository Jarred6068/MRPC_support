




pdf("/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/problem_trio_plots/WB_554_1407_plot.pdf")
A=ggplot(data=data.new, aes(x=cis.gene, y=trans.gene, color=SNP))+
  geom_point(size=2)+
  ggtitle(paste0(gn, collapse = ":"))
B=ggplot(data=data.new)+geom_boxplot(aes(x=SNP, y=cis.gene))
C=ggplot(data=data.new)+geom_boxplot(aes(x=SNP, y=trans.gene))
#D=ggplot(data=data.new)+geom_boxplot(aes(x=cis.gene, y=trans.gene))
sp=ggarrange(B, C, labels=c("B","C"), nrow=1, ncol=2)
sp2=ggarrange(A, sp, labels = c("A"), nrow = 2, ncol = 1)
plot(sp2)
dev.off()





pdf("/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/problem_trio_plots/WB_2789_7321_plot.pdf")
A=ggplot(data=data.new, aes(x=cis.gene, y=trans.gene, color=SNP))+
  geom_point(size=2)+
  ggtitle(paste0(gn, collapse = ":"))
B=ggplot(data=data.new)+geom_boxplot(aes(x=SNP, y=cis.gene))
C=ggplot(data=data.new)+geom_boxplot(aes(x=SNP, y=trans.gene))
#D=ggplot(data=data.new)+geom_boxplot(aes(x=cis.gene, y=trans.gene))
sp=ggarrange(B, C, labels=c("B","C"), nrow=1, ncol=2)
sp2=ggarrange(A, sp, labels = c("A"), nrow = 2, ncol = 1)
plot(sp2)
dev.off()


#read in postprocessing functions
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
#post process the significant mediation test pvalues by qvalue method



#cross analysis for GMAC inferred cis and trans trios (l1, l2, respectively) that did not match/overlap the 
#ADDIS inferred M1T1 and M1T2 trios
#run for WholeBlood
l1=cross.analyze(tissues=tissues.vec[,1], save=TRUE)
#view the breakdown of GMAC inferred significant trios across Tissues and Addis Model Class 
print(l1$breakdown.tab)

#look at first 10 entries of table for wholeblood
l1$final.tables[[5]][1:10,]


#Look at the addis inferred classes acr0ss at the intersection of the gmac inferred cis and trans models that 
#are unmatched with addis
ovlp.unm.gmac=l1$fl.unm[[1]][na.omit(row.match(l1$fl.unm[[1]], l2$fl.unm[[1]])), ]
summary(as.factor(ovlp.unm.gmac$Addis.infer.Class))/sum(summary(as.factor(ovlp.unm.gmac$Addis.infer.Class)))
#we notice that the majority are M3 and M4 types (fully connected and SNP regulation of cis/trans) ~80%

#create and save a histogram for the distribution of classes
png(filename = "/mnt/ceph/jarredk/GMACanalysis/GMAC.NO.ADDIS.WholeBlood.png")
plot.dist(class.vec.cis=l1$fl.unm[[1]][,5], class.vec.trans=l2$fl.unm[[1]][,5])
dev.off()


#checking for differences in numbers of selected PCs

#distribution of Number of PC's in matched trios
stem(rowSums(output.WB$cov.indicator.list[l1$fl.m[[1]]$index.in.trio.mat, ]))
#distribution of Number of PC's in unmatched trios
stem(rowSums(output.WB$cov.indicator.list[l1$fl.unm[[1]]$index.in.trio.mat, ]))
source("/mnt/ceph/jarredk/GMACanalysis/calc_meanpcs.R")
#summary of the mean number of pcs and sd across model types (cis/trans), for addis matching and non-matching inferred 
#trios and across tissues
print(mean.num.pc)
print(sd.num.pc)



#lookup the ADDIS and LOND inferred models for the first non-matching trio index
#outputs the correlation, lond, and addis edge matrices




print("#-------------------------------TRIO=9, Classed M4----------------------------------------")

l1$fl.unm[[1]][1,]

#Example 1 trio 9 class addis --> M4, class GMAC both M1T1 and M1T2
Lond2Addis.lookup(trio.index=9, tissue.name="WholeBlood", with.pc=TRUE)[2:4]

#look at the regressions for the two models:

list.data9=cross.regress(tissue="WholeBlood", trio.ind=9, mod.type="cis", addis.pcs=c("PC12","PC16"))

#From the regression on the trans gene following the GMAC process, we can see that both the SNP and cis gene are
#highly significant even in the presence of the adaptively selected confounders

#looking deeper: GMAC states that the null hypothesis assumes no correlation between the cis gene, Y, and the 
#trans gene, X, when conditioned on the eSNP, S, the alternative hypothesis is that mediation remains even
#after conditioning on S and permuting within each genotype. 

#condition the data on geno type and check the correlation between cis and trans genes within each group

#for cor(X,Y) | S = 0  --> homozygous reference

cor(list.data9$addis[which(list.data9$addis$SNP==0),][,1:2])

#for cor(X,Y) | S = 1  --> heterzygous reference

cor(list.data9$addis[which(list.data9$addis$SNP==1),][,1:2])

#for cor(X,Y) | S = 2  --> homozygous alternative

cor(list.data9$addis[which(list.data9$addis$SNP==2),][,1:2])

#We see here that there is a correlation between the cis and trans genes that is dependent on genotype
#e.g S=0, and S=2 (homozygous reference and alternative)
#ADDIS infers this model as fully connected meaning SNP regulates both genes and they are in turn associated
#with each other --> GMAC infers this model as mediation but does not directly account for a possible edge 
#between distinguish specifically between types of models. It determines if mediation exists between
#cis and trans in the presence of a SNP --> cis or SNP --> trans effect
fn="/mnt/ceph/jarredk/GMACanalysis/permreg_plots/GMAC.nomatch.permutedREG.trio9.png"
out.permreg=run.permuted.reg(trio=list.data9$GMAC, nperms=1000, plot=TRUE, filename=fn)

out.permreg$p.value
out.permreg$obs.wald.stat


r9=check.ns.gmac(input.list.data=list.data9, 
                 iters=1000, 
                 alpha=0.01, 
                 print.it=FALSE, 
                 which.mod="GMAC",
                 plot.it=TRUE, 
                 save.name="9M4")

r9a=check.ns.gmac(input.list.data=list.data9, 
                  iters=1000, 
                  alpha=0.01, 
                  print.it=FALSE, 
                  which.mod="ADDIS",
                  plot.it=TRUE, 
                  save.name="9M4")


print("#-------------------------------TRIO=11, Classed M4----------------------------------------")

l1$fl.unm[[1]][2,]

#Example 2 trio 11 class addis --> M4, class GMAC both M1T1 and M1T2
Lond2Addis.lookup(trio.index=11, tissue.name="WholeBlood", with.pc=TRUE)[2:4]

#look at the regressions for the two models:

list.data11=cross.regress(tissue="WholeBlood", trio.ind=11, mod.type="cis", addis.pcs=NULL)



#condition the data on geno type and check the correlation between cis and trans genes within each group

#for cor(X,Y) | S = 0  --> homozygous reference

cor(list.data11$addis[which(list.data11$addis$SNP==0),][,1:2])

#for cor(X,Y) | S = 1  --> heterzygous reference

cor(list.data11$addis[which(list.data11$addis$SNP==1),][,1:2])

#for cor(X,Y) | S = 2  --> homozygous alternative

cor(list.data11$addis[which(list.data11$addis$SNP==2),][,1:2])


fn="/mnt/ceph/jarredk/GMACanalysis/permreg_plots/GMAC.nomatch.permutedREG.trio11.png"
out.permreg=run.permuted.reg(trio=list.data11$GMAC, nperms=1000, plot=TRUE, filename=fn)

out.permreg$p.value
out.permreg$obs.wald.stat

r11=check.ns.gmac(input.list.data=list.data11, 
                  iters=1000, 
                  alpha=0.01, 
                  print.it=FALSE, 
                  which.mod="GMAC",
                  plot.it=TRUE, 
                  save.name="11M4")


r11a=check.ns.gmac(input.list.data=list.data11, 
                   iters=1000, 
                   alpha=0.01, 
                   print.it=FALSE, 
                   which.mod="ADDIS",
                   plot.it=TRUE, 
                   save.name="11M4")




print("#-------------------------------TRIO=13, Classed M4----------------------------------------")

l1$fl.unm[[1]][3,]

#Example 3 trio 13 class addis --> M4, class GMAC both M1T1 and M1T2
Lond2Addis.lookup(trio.index=13, tissue.name="WholeBlood", with.pc=TRUE)[2:4]

#look at the regressions for the two models:

list.data13=cross.regress(tissue="WholeBlood", trio.ind=13, mod.type="cis", addis.pcs=NULL)



#condition the data on geno type and check the correlation between cis and trans genes within each group

#for cor(X,Y) | S = 0  --> homozygous reference

cor(list.data13$addis[which(list.data13$addis$SNP==0),][,1:2])

#for cor(X,Y) | S = 1  --> heterzygous reference

cor(list.data13$addis[which(list.data13$addis$SNP==1),][,1:2])

#for cor(X,Y) | S = 2  --> homozygous alternative

cor(list.data13$addis[which(list.data13$addis$SNP==2),][,1:2])


fn="/mnt/ceph/jarredk/GMACanalysis/permreg_plots/GMAC.nomatch.permutedREG.trio13.png"
out.permreg=run.permuted.reg(trio=list.data13$GMAC, nperms=1000, plot=TRUE, filename=fn)

out.permreg$p.value
out.permreg$obs.wald.stat

r13=check.ns.gmac(input.list.data=list.data13, 
                  iters=1000, 
                  alpha=0.01, 
                  print.it=FALSE, 
                  which.mod="GMAC",
                  plot.it=TRUE, 
                  save.name="13M4")

r13a=check.ns.gmac(input.list.data=list.data13, 
                   iters=1000, 
                   alpha=0.01, 
                   print.it=FALSE, 
                   which.mod="ADDIS",
                   plot.it=TRUE, 
                   save.name="13M4")


print("#-------------------------------TRIO=17, Classed M3----------------------------------------")

l1$fl.unm[[1]][4,]

#Example 4 trio 17 class addis --> M3, class GMAC both M1T1 and M1T2
Lond2Addis.lookup(trio.index=17, tissue.name="WholeBlood", with.pc=TRUE)[2:4]

#look at the regressions for the two models:

list.data17=cross.regress(tissue="WholeBlood", trio.ind=17, mod.type="cis", addis.pcs=c("PC11","PC21","PC24"))



#condition the data on geno type and check the correlation between cis and trans genes within each group

#for cor(X,Y) | S = 0  --> homozygous reference

cor(list.data17$addis[which(list.data17$addis$SNP==0),][,1:2])

#for cor(X,Y) | S = 1  --> heterzygous reference

cor(list.data17$addis[which(list.data17$addis$SNP==1),][,1:2])

#for cor(X,Y) | S = 2  --> homozygous alternative

cor(list.data17$addis[which(list.data17$addis$SNP==2),][,1:2])


fn="/mnt/ceph/jarredk/GMACanalysis/permreg_plots/GMAC.nomatch.permutedREG.trio17.png"
out.permreg=run.permuted.reg(trio=list.data17$GMAC, nperms=1000, plot=TRUE, filename=fn)

out.permreg$p.value
out.permreg$obs.wald.stat


r17=check.ns.gmac(input.list.data=list.data17, 
                  iters=1000, 
                  alpha=0.01, 
                  print.it=FALSE, 
                  which.mod="GMAC",
                  plot.it=TRUE, 
                  save.name="17M3")

r17a=check.ns.gmac(input.list.data=list.data17, 
                   iters=1000, 
                   alpha=0.01, 
                   print.it=FALSE, 
                   which.mod="ADDIS",
                   plot.it=TRUE, 
                   save.name="17M3")





print("#-------------------------------TRIO=61, Classed M3----------------------------------------")

l1$fl.unm[[1]][10,]

#Example 4 trio 61 class addis --> M3, 
Lond2Addis.lookup(trio.index=61, tissue.name="WholeBlood", with.pc=TRUE)[2:4]

#look at the regressions for the two models:

list.data61=cross.regress(tissue="WholeBlood", trio.ind=61, mod.type="cis", addis.pcs=NULL)



#condition the data on geno type and check the correlation between cis and trans genes within each group

#for cor(X,Y) | S = 0  --> homozygous reference

cor(list.data61$addis[which(list.data61$addis$SNP==0),][,1:2])

#for cor(X,Y) | S = 1  --> heterzygous reference

cor(list.data61$addis[which(list.data61$addis$SNP==1),][,1:2])

#for cor(X,Y) | S = 2  --> homozygous alternative

cor(list.data61$addis[which(list.data61$addis$SNP==2),][,1:2])


fn="/mnt/ceph/jarredk/GMACanalysis/permreg_plots/GMAC.nomatch.permutedREG.trio17.png"
out.permreg=run.permuted.reg(trio=list.data61$GMAC, nperms=1000, plot=TRUE, filename=fn)

out.permreg$p.value
out.permreg$obs.wald.stat

#check numerical stability by 

r61=check.ns.gmac(input.list.data=list.data61, 
                  iters=1000, 
                  alpha=0.01, 
                  print.it=FALSE, 
                  which.mod="GMAC",
                  plot.it=TRUE, 
                  save.name="61M3")


r61a=check.ns.gmac(input.list.data=list.data61, 
                   iters=1000, 
                   alpha=0.01, 
                   print.it=FALSE, 
                   which.mod="ADDIS",
                   plot.it=TRUE, 
                   save.name="61M3")




print("#-------------------------------TRIO=12, Classed M1----------------------------------------")

l1$fl.m[[1]][1,]

#Example 1 of a matched trio: trio 12 class addis --> M1
Lond2Addis.lookup(trio.index=12, tissue.name="WholeBlood", with.pc=TRUE)[2:4]

#look at the regressions for the two models:

list.data12=cross.regress(tissue="WholeBlood", trio.ind=12, mod.type="cis", addis.pcs=NULL)



#condition the data on geno type and check the correlation between cis and trans genes within each group

#for cor(X,Y) | S = 0  --> homozygous reference

cor(list.data12$addis[which(list.data12$addis$SNP==0),][,1:2])

#for cor(X,Y) | S = 1  --> heterzygous reference

cor(list.data12$addis[which(list.data12$addis$SNP==1),][,1:2])

#for cor(X,Y) | S = 2  --> homozygous alternative

cor(list.data12$addis[which(list.data12$addis$SNP==2),][,1:2])


fn="/mnt/ceph/jarredk/GMACanalysis/permreg_plots/GMAC.nomatch.permutedREG.trio11.png"
out.permreg=run.permuted.reg(trio=list.data12$GMAC, nperms=1000, plot=TRUE, filename=fn)

out.permreg$p.value
out.permreg$obs.wald.stat

r12=check.ns.gmac(input.list.data=list.data12, 
                  iters=1000, 
                  alpha=0.01, 
                  print.it=FALSE, 
                  which.mod="GMAC",
                  plot.it=TRUE, 
                  save.name="12M1")


r12a=check.ns.gmac(input.list.data=list.data12, 
                   iters=1000, 
                   alpha=0.01, 
                   print.it=FALSE, 
                   which.mod="ADDIS",
                   plot.it=TRUE, 
                   save.name="12M1")





print("#-------------------------------TRIO=23, Classed M1----------------------------------------")

l1$fl.m[[1]][2,]

#Example 2 of a matched trio: trio 23 class addis --> M1
Lond2Addis.lookup(trio.index=23, tissue.name="WholeBlood", with.pc=TRUE)[2:4]

#look at the regressions for the two models:

list.data23=cross.regress(tissue="WholeBlood", trio.ind=23, mod.type="cis", addis.pcs=NULL)



#condition the data on geno type and check the correlation between cis and trans genes within each group

#for cor(X,Y) | S = 0  --> homozygous reference

cor(list.data23$addis[which(list.data23$addis$SNP==0),][,1:2])

#for cor(X,Y) | S = 1  --> heterzygous reference

cor(list.data23$addis[which(list.data23$addis$SNP==1),][,1:2])

#for cor(X,Y) | S = 2  --> homozygous alternative

cor(list.data23$addis[which(list.data23$addis$SNP==2),][,1:2])


fn="/mnt/ceph/jarredk/GMACanalysis/permreg_plots/GMAC.nomatch.permutedREG.trio11.png"
out.permreg=run.permuted.reg(trio=list.data23$GMAC, nperms=1000, plot=TRUE, filename=fn)

out.permreg$p.value
out.permreg$obs.wald.stat

r23=check.ns.gmac(input.list.data=list.data23, 
                  iters=1000, 
                  alpha=0.01, 
                  print.it=FALSE, 
                  which.mod="GMAC",
                  plot.it=TRUE, 
                  save.name="23M1")


r23a=check.ns.gmac(input.list.data=list.data23, 
                   iters=1000, 
                   alpha=0.01, 
                   print.it=FALSE, 
                   which.mod="ADDIS",
                   plot.it=TRUE, 
                   save.name="23M1")





print("#-------------------------------TRIO=87, Classed M1----------------------------------------")

l1$fl.m[[1]][4,]

#Example 3 of a matched trio: trio 87 class addis --> M1
Lond2Addis.lookup(trio.index=87, tissue.name="WholeBlood", with.pc=TRUE)[2:4]

#look at the regressions for the two models:

list.data87=cross.regress(tissue="WholeBlood", trio.ind=87, mod.type="cis", addis.pcs=NULL)



#condition the data on geno type and check the correlation between cis and trans genes within each group

#for cor(X,Y) | S = 0  --> homozygous reference

cor(list.data87$addis[which(list.data87$addis$SNP==0),][,1:2])

#for cor(X,Y) | S = 1  --> heterzygous reference

cor(list.data87$addis[which(list.data87$addis$SNP==1),][,1:2])

#for cor(X,Y) | S = 2  --> homozygous alternative

cor(list.data87$addis[which(list.data87$addis$SNP==2),][,1:2])


fn="/mnt/ceph/jarredk/GMACanalysis/permreg_plots/GMAC.nomatch.permutedREG.trio11.png"
out.permreg=run.permuted.reg(trio=list.data87$GMAC, nperms=1000, plot=TRUE, filename=fn)

out.permreg$p.value
out.permreg$obs.wald.stat

r87=check.ns.gmac(input.list.data=list.data87, 
                  iters=1000, 
                  alpha=0.01, 
                  print.it=FALSE, 
                  which.mod="GMAC",
                  plot.it=TRUE, 
                  save.name="23M1")


r87a=check.ns.gmac(input.list.data=list.data87, 
                   iters=1000, 
                   alpha=0.01, 
                   print.it=FALSE, 
                   which.mod="ADDIS",
                   plot.it=TRUE, 
                   save.name="23M1")


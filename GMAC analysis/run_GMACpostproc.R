
#read in postprocessing functions
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")

#this function compares inferred trios in GMAC and 
out1=match.trios(tissues=tissues.vec[,1], which.mrpc="ADDIS")
#this table shows the overlapping cis and trans mediation trios between ADDIS and GMAC
out1$table
#the first five non-overlapping trios inferred to be cis mediation by GMAC
out1$unmatched.cis[[1]][1:5,]
#print dimension
print(dim(out1$unmatched.cis[[1]]))
#the first five non-overlapping trios inferred to be trans mediation by GMAC
out1$unmatched.trans[[1]][1:5,]
#print dimension
print(dim(out1$unmatched.trans[[1]]))
#count the intersect of the unmatched cis and trans trios
row.match(out1$unmatched.cis[[1]],out1$unmatched.trans[[1]])

#run now for LOND
# out2=match.trios(tissues=tissues.vec[,1], which.mrpc="LOND")
# #this table shows the overlapping cis and trans mediation trios between LOND and GMAC
# out2$table
# #the first five non-overlapping trios inferred to be cis mediation by GMAC
# out2$unmatched.cis[[1]][1:5,]
# #the first five non-overlapping trios inferred to be trans mediation by GMAC
# out2$unmatched.trans[[1]][1:5,]


#cross analysis for GMAC inferred cis and trans trios (l1, l2, respectively) that did not match/overlap the 
#ADDIS inferred M1T1 and M1T2 trios
#run for WholeBlood
l1=cross.analyze(tissues="WholeBlood", which.type="cis", save=TRUE)
#view the first 10 non-overlapping cis mediated trios and their model type under ADDIS
l1[[1]][1:10,]

l2=cross.analyze(tissues="WholeBlood", which.type="trans", save=TRUE)
#view the first 10 non-overlapping trans trios and their model type under ADDIS
l2[[1]][1:10,]

#what is the distribution of addis model types that overlap between both GMAC inferred cis/trans mediation 
#trios
ovlp.gmac=l1[[1]][na.omit(row.match(l1[[1]], l2[[1]])), ]
summary(as.factor(ovlp.gmac$Addis.infer.Class))/sum(summary(as.factor(ovlp.gmac$Addis.infer.Class)))
#we notice that the majority are M3 and M4 types (fully connected and SNP regulation of cis/trans)

#create and save a histogram for the distribution of classes
png(filename = "/mnt/ceph/jarredk/GMACanalysis/GMAC.NO.ADDIS.WholeBlood.png")
plot.dist(class.vec.cis=l1[[1]][,5], class.vec.trans=l2[[1]][,5])
dev.off()

#lookup the ADDIS and LOND inferred models for the first non-matching trio index
#outputs the correlation, lond, and addis edge matrices




#-------------------------------TRIO=9, Classed M4----------------------------------------

l1[[1]][1,]

#Example 1 trio 9 class addis --> M4, class GMAC both M1T1 and M1T2
Lond2Addis.lookup(trio.index=9, tissue.name="WholeBlood", with.pc=TRUE)[2:4]

#look at the regressions for the two models:

list.data=cross.regress(tissue="WholeBlood", trio.ind=9, mod.type="cis", addis.pcs=c("PC12","PC16"))

#From the regression on the trans gene following the GMAC process, we can see that both the SNP and cis gene are
#highly significant even in the presence of the adaptively selected confounders

#looking deeper: GMAC states that the null hypothesis assumes no correlation between the cis gene, Y, and the 
#trans gene, X, when conditioned on the eSNP, S, the alternative hypothesis is that mediation remains even
#after conditioning on S and permuting within each genotype. 

#condition the data on geno type and check the correlation between cis and trans genes within each group

#for cor(X,Y) | S = 0  --> homozygous reference

cor(list.data$addis[which(list.data$addis$SNP==0),][,1:2])

#for cor(X,Y) | S = 1  --> heterzygous reference

cor(list.data$addis[which(list.data$addis$SNP==1),][,1:2])

#for cor(X,Y) | S = 2  --> homozygous alternative

cor(list.data$addis[which(list.data$addis$SNP==2),][,1:2])

#We see here that there is a correlation between the cis and trans genes that is dependent on genotype
#e.g S=0, and S=2 (homozygous reference and alternative)
#ADDIS infers this model as fully connected meaning SNP regulates both genes and they are in turn associated
#with each other --> GMAC infers this model as mediation but does not directly account for a possible edge 
#between the trans gene and the SNP. 
fn="/mnt/ceph/jarredk/GMACanalysis/permreg_plots/GMAC.nomatch.permutedREG.trio9.png"
out.permreg=run.permuted.reg(trio=list.data$GMAC, nperms=10000, plot=TRUE, filename=fn)

out.permreg$p.value
out.permreg$obs.wald.stat


#-------------------------------TRIO=11, Classed M4----------------------------------------

l1[[1]][2,]

#Example 1 trio 9 class addis --> M4, class GMAC both M1T1 and M1T2
Lond2Addis.lookup(trio.index=11, tissue.name="WholeBlood", with.pc=TRUE)[2:4]

#look at the regressions for the two models:

list.data=cross.regress(tissue="WholeBlood", trio.ind=11, mod.type="cis", addis.pcs=NULL)



#condition the data on geno type and check the correlation between cis and trans genes within each group

#for cor(X,Y) | S = 0  --> homozygous reference

cor(list.data$addis[which(list.data$addis$SNP==0),][,1:2])

#for cor(X,Y) | S = 1  --> heterzygous reference

cor(list.data$addis[which(list.data$addis$SNP==1),][,1:2])

#for cor(X,Y) | S = 2  --> homozygous alternative

cor(list.data$addis[which(list.data$addis$SNP==2),][,1:2])


fn="/mnt/ceph/jarredk/GMACanalysis/permreg_plots/GMAC.nomatch.permutedREG.trio11.png"
out.permreg=run.permuted.reg(trio=list.data$GMAC, nperms=10000, plot=TRUE, filename=fn)

out.permreg$p.value
out.permreg$obs.wald.stat






#-------------------------------TRIO=13, Classed M4----------------------------------------

l1[[1]][3,]

#Example 1 trio 9 class addis --> M4, class GMAC both M1T1 and M1T2
Lond2Addis.lookup(trio.index=13, tissue.name="WholeBlood", with.pc=TRUE)[2:4]

#look at the regressions for the two models:

list.data=cross.regress(tissue="WholeBlood", trio.ind=13, mod.type="cis", addis.pcs=NULL)



#condition the data on geno type and check the correlation between cis and trans genes within each group

#for cor(X,Y) | S = 0  --> homozygous reference

cor(list.data$addis[which(list.data$addis$SNP==0),][,1:2])

#for cor(X,Y) | S = 1  --> heterzygous reference

cor(list.data$addis[which(list.data$addis$SNP==1),][,1:2])

#for cor(X,Y) | S = 2  --> homozygous alternative

cor(list.data$addis[which(list.data$addis$SNP==2),][,1:2])


fn="/mnt/ceph/jarredk/GMACanalysis/permreg_plots/GMAC.nomatch.permutedREG.trio11.png"
out.permreg=run.permuted.reg(trio=list.data$GMAC, nperms=10000, plot=TRUE, filename=fn)

out.permreg$p.value
out.permreg$obs.wald.stat






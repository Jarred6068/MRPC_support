
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


#create and save a histogram for the distribution of classes
png(filename = "/mnt/ceph/jarredk/GMACanalysis/GMAC.NO.ADDIS.WholeBlood.png")
plot.dist(class.vec.cis=l1[[1]][,5], class.vec.trans=l2[[1]][,5])
dev.off()

#lookup the ADDIS and LOND inferred models for the first non-matching trio index
#outputs the correlation, lond, and addis edge matrices
Lond2Addis.lookup(trio.index=9, tissue.name="WholeBlood", with.pc=FALSE)[2:4]

#look at the regressions for the two models:

cross.regress(tissue="WholeBlood", trio.ind=9, mod.type="cis", addis.pcs=NULL)





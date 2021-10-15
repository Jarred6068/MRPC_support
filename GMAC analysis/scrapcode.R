
  #read in postprocessing functions

source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
l1=cross.analyze(tissues=tissues.vec[5,1], save=FALSE)
#check first 15 M3's for wholeblood
data1=l1$final.tables[[1]][which(l1$final.tables[[1]]$Addis.Class=="M3"),]
M3.trios.sample=l1$final.tables[[1]][which(l1$final.tables[[1]]$Addis.Class=="M3"),]$Trio.Num

#Lond2Addis.lookup(trio.index=135, tissue.name="WholeBlood", with.pc=TRUE)[2:4]
outM3=runit(indata = data1[3,], trio.number =  M3.trios.sample[3], mtype = "M3")



test.data=outM3$datalist$GMAC[,-c((dim(outM3$datalist$GMAC)[2]-2):dim(outM3$datalist$GMAC)[2])]
test.data=cbind.data.frame(SNP=test.data[,3],test.data[,1:2], test.data[,4:dim(test.data)[2]])

n <- nrow (test.data)
V <- colnames(test.data)     # Column names

# Classical correlation
suffStat <- list(C = cor(test.data, use = "complete.obs"),
                 n = n)
#run MRPC on TRIO with pc's
MRPC.fit.FDR.addis <- MRPC(test.data,
                           suffStat,
                           GV = 1,
                           FDR = 0.05,
                           indepTest = 'gaussCItest',
                           labels = V,
                           verbose = FALSE)

as(MRPC.fit.FDR.addis@graph, "matrix")[1:3,1:3]


Lond2Addis.lookup(trio.index=17, tissue.name="WholeBlood", with.pc=TRUE)

list.data=cross.regress(tissue="WholeBlood", trio.ind=17, mod.type="cis", addis.pcs=NULL)


fname1=paste0("/mnt/ceph/jarredk/GMACanalysis/additional_plots/GMAC.nomatch.permutedREG.trio",
              319,".png")
out.permreg=run.permuted.reg(trio=list.data$GMAC, nperms=10000, plot=TRUE, filename=fname1)


#========================================================================================
##suppression analysis
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
l1=cross.analyze(tissues=tissues.vec[5,1], save=FALSE)
all.m3=l1$final.tables[[1]][which(l1$final.tables[[1]]$Addis.Class=="M3"),]
print(all.m3)
list.data=cross.regress(tissue="WholeBlood", trio.ind=17, mod.type="cis", addis.pcs=NULL)
#get data & convert factors
test.data=list.data$GMAC
test.data$pcr=as.factor(test.data$pcr)
test.data$platform=as.factor(test.data$platform)
test.data$sex=as.factor(test.data$sex)

test.mod=lm(trans.gene~., data=test.data)
vif(test.mod)

#look at suppression index for the first 2 independent variables
id.suppress(list.data$GMAC[,1:7], verbose=TRUE)












nn=10
l1=cross.analyze(tissues="WholeBlood", save=FALSE)
all.m0m3=l1$final.tables[[1]][c(which(l1$final.tables[[1]]$Addis.Class=="M3"), which(l1$final.tables[[1]]$Addis.Class=="MO")),]
trios=sample(all.m0m3$Trio.Num, nn)
ot1=run.simu12(tissue = "WholeBlood" ,trios=trios, 
               mod.type.vec=all.m0m3$Mediation.type[match(trios, all.m0m3$Trio.Num)],
               alpha=0.001, n=10)



list.data=cross.regress(tissue="WholeBlood", trio.ind=121, mod.type="trans", addis.pcs=NULL)
simu2(tissue = "WholeBlood", data=list.data$GMAC, seed=NULL, mod.type="Trans.Med", n=10, verbose=TRUE)
simu1(data=list.data$GMAC,alpha=0.001, mod.type=NULL, verbose=TRUE)

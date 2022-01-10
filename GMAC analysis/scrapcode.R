
  #read in postprocessing functions

source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
l1=cross.analyze(tissues="CellsEBVtransformedlymphocytes", save=FALSE)
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

Lond2Addis.lookup(trio.index=54, tissue.name="AdiposeSubcutaneous", with.pc=TRUE)

list.data=cross.regress(tissue="AdiposeSubcutaneous", trio.ind=54, mod.type="both", addis.pcs=NULL)

run.permuted.reg(trio=list.data$GMAC, nperms=1000, plot=FALSE, filename=NULL, Alg="GMAC", med.type="Trans.Med")

fname1=paste0("/mnt/ceph/jarredk/GMACanalysis/GMAC.permutedREG.trio", 8112,".png")
start.t=Sys.time()
out.permreg=run.permuted.reg(trio=list.data$GMAC, nperms=1000, plot=F, filename=fname1)
end.t=Sys.time()
end.t-start.t


perm6=run.permuted.reg(sim$sim.data, nperms=1000, Alg="GMAC", med.type = "Both")

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
ot1=run.simu12(tissue = "WholeBlood" ,trios=trios, l1.table=all.m0m3,
               mod.type.vec=all.m0m3$Mediation.type[match(trios, all.m0m3$Trio.Num)],
               alpha=0.001, n="random", seed=NULL)

run.permuted.reg(list.data$addis, nperms=1000, plot=FALSE, filename=NULL, Alg="ADDIS")$p.value

list.data=cross.regress(tissue="AdiposeSubcutaneous", trio.ind=5683, mod.type="cis", addis.pcs=NULL)
x=simu2(tissue = "WholeBlood", data=list.data$GMAC, seed=NULL, mod.type="Both", n=10, verbose=TRUE)
xx=simu1(data=list.data$GMAC,alpha=0.001, mod.type="Both", verbose=TRUE)

x=simu2(tissue = "WholeBlood", data=list.data$GMAC, seed=NULL, mod.type="Both", n="random", verbose=TRUE)

#generate table
for(i in 1:length(tissues.vec[,1])){
  l1=cross.analyze(tissues=tissues, save=FALSE)
  
}




run.permuted.reg(trio=list.data$GMAC, nperms=1000, plot=TRUE, filename="/mnt/ceph/jarredk/GMACanalysis/trio5683ASpermuted.png", Alg="GMAC")


check.ns.gmac(input.list.data=list.data, iters=1, alpha=0.01, print.it=FALSE, which.mod="GMAC",
                       plot.it=T, save.name="TEST5683AS_CHECK1")


#looking at methyl cors
V=NULL
for(i in c(1:22, "X", "Y")){
  A=read.table(file=paste0("/mnt/ceph/jarredk/Methyl/BSGS_Correlations/BSGS_Correlations_20210914/M_chr", 
                           i,"_correlations.txt"), sep="\t", header=T)
  #print(paste0("chr",i,"=",length(unique(A$Methyl_Probe_ID[which(is.na(A$Cor))]))))
  print(paste0("Unique Probes for chr",i,"=",length(unique(A$Methyl_Probe_ID))))
  
  V[i]=length(unique(A$Methyl_Probe_ID))
  
}

print(sum(V))








#--------------------------------------------------------------------------------------------------

#Add in simulation for MRPC on GMAC model to the simulations.txt dataset
#source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
sim=read.table(file=paste0("/mnt/ceph/jarredk/GMACanalysis/master_tables/simulations_","AdiposeSubcutaneous",".txt"), sep="\t", header = T)

p=NULL

for(i in 1:dim(sim)[1]){
  print(i)
  
  L2A=Lond2Addis.lookup(trio.index=sim$Trio.Num[i], tissue.name=sim$Tissue[i], with.pc=TRUE)
  
  if(length(colnames(L2A$correlation))>3){
    addis.pcs=colnames(L2A$correlation)[-c(1:3)]
  }else{
    addis.pcs=NULL
  }
  
  list.data=cross.regress(tissue=sim$Tissue[i], trio.ind=sim$Trio.Num[i], mod.type=sim$Med.type[i], addis.pcs=addis.pcs, verbose=F)
  
  out3=simu3(MRPC.data=list.data$addis, GMAC.data=list.data$GMAC, mod.type=sim$Med.type[i], verbose=F)
  
  p[i]=out3$pvalue
  
}

sim$TGM.p=p












#checking trios with conflicting mediation and permutation p values
#must be run post run_simulations.slurm on cluster
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
sim=read.table(file="/mnt/ceph/jarredk/GMACanalysis/master_tables/simulations_adipose_subcutaneous.txt", sep="\t", header = T)

#wt=sim[which(sim$Med.pvalue.GMAC>0.05),]

snp.dist=as.data.frame(matrix(0, nrow = dim(sim)[1], ncol = 3))
colnames(snp.dist)=c("hr","het","ha")

for(i in 1:dim(sim)[1]){
  
  if(sim$Med.type[i]=="Trans.Med"){
    list.data=cross.regress(tissue=sim$Tissue[i], trio.ind=sim$Trio.Num[i], mod.type="trans", addis.pcs=NULL, verbose=F)
  }else{
    list.data=cross.regress(tissue=sim$Tissue[i], trio.ind=sim$Trio.Num[i], mod.type="cis", addis.pcs=NULL, verbose=F)
  }
  
  
  #account for SNPs with less than 3 genotypes
  if(length(summary(factor(list.data$GMAC$SNP)))<3){
    
    if( length(which(names(summary(factor(list.data$GMAC$SNP)))=="0"))==0 ){
      snp.dist[i,]=c(NA ,summary(factor(list.data$GMAC$SNP)))
    }else if ( length(which(names(summary(factor(list.data$GMAC$SNP)))=="1"))==0 ){
      snp.dist[i,]=c(summary(factor(list.data$GMAC$SNP))[1], NA , summary(factor(list.data$GMAC$SNP))[2])
    }else{
      snp.dist[i,]=c(summary(factor(list.data$GMAC$SNP)), NA)
    }
    
  }else{
    
    snp.dist[i,]=summary(factor(list.data$GMAC$SNP))
    
  }
  
}


output1=cbind.data.frame(sim, snp.dist)

write.csv(output1, file = "/mnt/ceph/jarredk/GMACanalysis/master_tables/TRIOS_imbalanced_genotypes_AS.csv", row.names = F)







#Add in simulation for MRPC on GMAC model to the simulations.txt dataset
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
sim=read.table(file="/mnt/ceph/jarredk/GMACanalysis/master_tables/simulations_adipose_subcutaneous.txt", sep="\t", header = T)

p=NULL

for(i in 1:dim(sim)[1]){
  print(i)
  
  L2A=Lond2Addis.lookup(trio.index=sim$Trio.Num[i], tissue.name=sim$Tissue[i], with.pc=TRUE)
  
  if(length(colnames(L2A$correlation))>3){
    addis.pcs=colnames(L2A$correlation)[-c(1:3)]
  }else{
    addis.pcs=NULL
  }
  
  list.data=cross.regress(tissue=sim$Tissue[i], trio.ind=sim$Trio.Num[i], mod.type=sim$Med.type[i], addis.pcs=addis.pcs, verbose=F)
  
  p[i]=simu3(MRPC.data=list.data$addis, GMAC.data=list.data$GMAC, mod.type=sim$Med.type[i], verbose=F)$pvalue
  
  
}

sim$TGM.p=p

write.table(sim, file = "/mnt/ceph/jarredk/GMACanalysis/master_tables/simulations_adipose_subcutaneous.txt", col.names = T, row.names = F, sep = "\t",
            quote=F)






#Add in permutation pvalue for simulation for MRPC on GMAC model to the simulations.txt dataset
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
t1=read.csv(file = "/mnt/ceph/jarredk/GMACanalysis/master_tables/TRIOS_imbalanced_genotypes_WB.csv")

p=NULL


for(i in 1:dim(t1)[1]){
  #print(i)
  
  L2A=Lond2Addis.lookup(trio.index=t1$Trio.Num[i], tissue.name=t1$Tissue[i], with.pc=TRUE)
  print(dim(L2A$correlation))
  if(length(colnames(L2A$correlation))>3){
    addis.pcs=colnames(L2A$correlation)[-c(1:3)]
  }else{
    addis.pcs=NULL
  }
  
  list.data=cross.regress(tissue=t1$Tissue[i], trio.ind=t1$Trio.Num[i], mod.type=t1$Med.type[i], addis.pcs=addis.pcs, verbose=F)
  
  out=simu3(MRPC.data=list.data$addis, GMAC.data=list.data$GMAC, mod.type=t1$Med.type[i], verbose=F)
  print(dim(out$sim.data))
  p[i]=run.permuted.reg(out$sim.data, nperms=1000, Alg="ADDIS", med.type = t1$Med.type[i])$p.value
  print(p[i])
  
}


t1$TGM.perm=p


write.csv(t1, file = "/mnt/ceph/jarredk/GMACanalysis/master_tables/TRIOS_imbalanced_genotypes_WB.csv", row.names = F)



























#testing GMAC inference on trios with rare alleles in cis-eQTL:
source("/mnt/ceph/jarredk/GMACanalysis/GMACanalysis.R")
library('GMAC', lib='/mnt/ceph/jarredk/Rpackages')
library('mice', lib='/mnt/ceph/jarredk/Rpackages')
library('missMDA', lib='/mnt/ceph/jarredk/Rpackages')
sim.odd=read.csv(file="mnt/ceph/jarredk/GMACanalysis/master_tables/TRIOS_imbalanced_genotypes.csv", header = T)
#trios.to.keep=
t=1

#load additional files
file1=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",
            tissues.vec[t,1],"_AllPC/data.snp.cis.trans.final.",
            tissues.vec[t,1],".V8.unique.snps.RData", sep = "")

file2=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",
            tissues.vec[t,1],"_AllPC/PCs.matrix.",
            tissues.vec[t,1], ".RData", sep="")

#load trio, pc, and confounders data
trios=loadRData(fileName=file1)
PCs=loadRData(fileName=file2)
#sig.asso.pcs=loadRData(fileName=file3)
#edata=loadRData(fileName=file4)
confounders=read.table(paste0("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/GTEx_Analysis_v8_eQTL_covariates/", 
                              tissues.vec[t,2], '.v8.covariates.txt'), header = T, sep="\t", row.names = 1)[66:68,]

print(paste('data loaded for', tissues.vec[t,2], sep = ' '))

#run assemble tables function
tables.gmac.list=assemble.tables(trio.table=trios, PCmat = PCs, kc = confounders, which.imp = 'miss')

output <- gmac(known.conf = tables.gmac.list$known.conf, cov.pool = tables.gmac.list$cov.pool, 
               exp.dat = tables.gmac.list$exp.dat[row.match(sim.odd[1,18:19], matrix(row.names(tables.gmac.list$exp.dat)), ncol=2, byrow=T ),],
               snp.dat.cis = tables.gmac.list$snp.dat.cis[5863,], 
               trios.idx = tables.gmac.list$trios.idx[5683,], nperm = 10000, nominal.p = TRUE)


















#---------------------simulation-tables:-------------------------

source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")


#nn=50

#small truth model large inferred model




l1=cross.analyze(tissues="WholeBlood", save=FALSE)

trios=l1$final.tables[[1]]$Trio.Num
med.types=l1$final.tables[[1]]$Mediation.type
print("running Simulations...")

start.time=Sys.time()
ot1=run.simu12(tissue = "WholeBlood" ,trios=trios[1:10], 
               l1.table=l1$final.tables[[1]][1:10,],
               mod.type.vec=med.types[1:10],
               alpha=0.001, n="random")
end.time=Sys.time()
print("...done")
print(paste0("time = ", end.time-start.time))
ot1$Tissue=rep("WholeBlood")





















l1=cross.analyze(tissues="AdiposeSubcutaneous", save=FALSE)

source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")

list.data=cross.regress(tissue="AdiposeSubcutaneous", trio.ind=9711, mod.type="trans", addis.pcs=NULL)

run.permuted.reg(trio=list.data$GMAC, nperms=1000, plot=FALSE, filename=NULL, Alg="GMAC", med.type="Trans.Med")$null.wald.stat
















#code to run multiple simulations to check the stability of the analysis

#load in needed stuff
sim=read.csv(file="/mnt/ceph/jarredk/GMACanalysis/master_tables/TRIOS_imbalanced_genotypes_WB.csv", header = T)
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")

tissue="WholeBlood"
l1=cross.analyze(tissues=tissue, save=FALSE)
#pull out trios which were no longer significant after the GSS simulation 
#trios=sim$Trio.Num[sim$Nominal.p.GSS>0.05][1:5]
#mod.type.vec=sim$Med.type[sim$Nominal.p.GSS>0.05][1:5]
trios = sim$Trio.Num
mod.type.vec = sim$Med.type

ss=1000
storage1=vector("list", length = length(trios))
plot.it=FALSE
STM.median.p=NULL
LTM.median.p=NULL
TGM.median.p=NULL
GSS.median.p=NULL

for(i in 1:5){
  
  
  print(paste0("getting data for trio # ", trios[i]))
  if(mod.type.vec[i]=="Both"){
    
    list.data=cross.regress(tissue=tissue, 
                            trio.ind=356, 
                            mod.type="cis", 
                            #addis.pcs=addis.pcs, 
                            verbose = TRUE)
    
  }else if(mod.type.vec[i]=="Cis.Med"){
    
    list.data=cross.regress(tissue=tissue, 
                            trio.ind=trios[i], 
                            mod.type="cis", 
                            #addis.pcs=addis.pcs, 
                            verbose = FALSE)
  }else{
    
    list.data=cross.regress(tissue=tissue, 
                            trio.ind=trios[i], 
                            mod.type="trans", 
                            #addis.pcs=addis.pcs, 
                            verbose = FALSE)
    
  }
  print("...done...")
  
  p1=as.data.frame(matrix(0, nrow = ss, ncol = 4))
  b1=as.data.frame(matrix(0, nrow = ss, ncol = 4))
  colnames(p1)=c("simu1.p", "simu2.p", "simu3.p", "simu4.p")
  colnames(b1)=c("simu1.b", "simu2.b", "simu3.b", "simu4.b")
  
  print(paste0("running all simulations for trio # ", trios[i]))
  for(j in 1:ss){
    
    
    out=simu1(list.data$GMAC, alpha = 0.001, mod.type=mod.type.vec[i], verbose=F)
    nn=floor(runif(1, 1, 20))
    #print(nn)
    out2=simu2(tissue = tissue, data=list.data$GMAC, mod.type=mod.type.vec[i], n=nn, verbose=F)
    out3=simu3(MRPC.data=list.data$addis, GMAC.data=list.data$GMAC, mod.type=mod.type.vec[i], verbose=F)
    out4=simu4(GMAC.data=list.data$GMAC, mod.type=mod.type.vec[i], verbose=F)

    p1[j,1] = out$pvalue
    p1[j,2] = out2$pvalue
    p1[j,3] = out3$pvalue
    p1[j,4] = out4$pvalue
    
    b1[j,1] = out$b.coef
    b1[j,2] = out2$b.coef
    b1[j,3] = out3$b.coef
    b1[j,4] = out4$b.coef

  }
  print("...done!")
  
  STM.median.p[i] = median(p1$simu1.p)
  LTM.median.p[i] = median(p1$simu2.p)
  TGM.median.p[i] = median(p1$simu3.p)
  GSS.median.p[i] = median(p1$simu4.p)
  
  
  storage1[[i]] = cbind.data.frame(p1, b1)
  
  
  if(plot.it==TRUE){
    png(paste0("/mnt/ceph/jarredk/GMACanalysis/Unstable_trios_plots/plot_", tissue, "_", trios[i], "_",mod.type.vec[i], ".png"))
    plot(p1, p3, pch=21, bg="black",
         xlab="Simulated Para. P-value",
         ylab="Simulated Nominal P-value",
         main = paste0("GSS 1000 simulations of trio = ", trios[i], ";", tissue ))
    abline(a=0, b=1, lty="dotted", col="red")
    dev.off()
  
    png(paste0("/mnt/ceph/jarredk/GMACanalysis/Unstable_trios_plots/histogram_", tissue, "_", trios[i], "_",mod.type.vec[i], ".png"))
    par(mfrow=c(2,2))
    hist(p1, xlab="Simulated Para. P-value",
         ylab="Simulated Nominal P-value",
         main = paste0("GSS 1000 simulations of trio = ", trios[i], ";", tissue ))
    hist(p2,
         xlab="Simulated Para. P-value",
         ylab="Simulated Nominal P-value",
         main = paste0("GSS 1000 simulations of trio = ", trios[i], ";", tissue ))
    hist(p3,
         xlab="Simulated Para. P-value",
         ylab="Simulated Nominal P-value",
         main = paste0("GSS 1000 simulations of trio = ", trios[i], ";", tissue ))
    dev.off()
  }
}


sim.new=cbind.data.frame(sim, STM.median.p, LTM.median.p, TGM.median.p, GSS.median.p)

save(storage1, file = "/mnt/ceph/jarredk/GMACanalysis/master_tables/multi_sim_data.RData")
write.csv(sim.new, file = "/mnt/ceph/jarredk/GMACanalysis/master_tables/TRIOS_imbalanced_genotypes_WB_Updated_12_10_2021.csv")











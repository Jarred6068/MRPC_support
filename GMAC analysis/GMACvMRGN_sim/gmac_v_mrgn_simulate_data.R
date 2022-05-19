

#Simulate All Data and Parameters:
library("MRGN",lib="/mnt/ceph/jarredk/Rpackages")
library("GMAC", lib="/mnt/ceph/jarredk/Rpackages")
source("/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/simulate_trios_with_c_kc.R")
#set.seed(234)
#pre-allocate simulation conditions
model.types=c("model0","model1","model2","model3","model4")
number.of.datasets=10000
#preallocate
params=as.data.frame(matrix(0, nrow = number.of.datasets, ncol = 6))
colnames(params)=c("model","SD", "minor.freq", "b.snp","b.med","number.of.conf")
#load in pc mat from whole blood
tissue.name="WholeBlood"
pc.matrix=loadRData(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissue.name,
                          "_AllPC/PCs.matrix.",tissue.name,".RData", sep = ""))
#list for all simulated datasets
small.datasets=list()


kc = as.matrix(cbind(pcr = sample(c(1:3), size = dim(pc.matrix)[1], replace = T, prob = rep(1/3,3)),
           sex = rbinom(dim(pc.matrix)[1], 1, 0.5)+1,
           age = rnorm(dim(pc.matrix)[1], mean = 35, sd=5)))


print("Simulating Data...")
#pc.matrix=as.matrix(pc.matrix)
#simulate all datasets
resamp.vec=NULL
for(j in 1:number.of.datasets){
  init=0
  k=1
  resamples=NULL
  #simulate trio under parameters
  while(length(init)<=2){

    #simulate parameters
    params$model[j]=sample(model.types, 1)
    params$minor.freq[j]=sample(seq(0.01, 0.5, 0.01),1)
    params$b.snp[j]=sample(seq(0.5, 1.5, 0.1), 1)
    params$b.med[j]=sample(seq(0.5, 1, 0.1), 1)
    params$SD[j]=(sample(seq(0.3, 1.5, 0.1), 1)*params$b.med[j])
    u = abs(round(rnorm(1, 10, 2)))
    while(u<=0){u = abs(round(rnorm(1, 10, 2)))} #so that average is 10 and min = 2 to %*% error
    #simulate data under parameters
    X=simData(theta=params$minor.freq[j],
              model=params$model[j],
              b0.1=0,
              b.snp=params$b.snp[j],
              b.med=params$b.med[j],
              sd.1=params$SD[j],
              G = pc.matrix,
              u = u,
              kc = kc)
    params$number.of.conf[j]=u
    #params$max.conf.sig[j]= w
    init=unique(X$V1)
    resamples[k]=k
    k=k+1
  }

  #print(paste0(length(resamples)," before all 3 genotypes represented"))
  #print(dim(X))
  #store
  print(params[j,])
  print(dim(X))
  small.datasets[[j]]=X
}

print("Done!...Saving...")

save(small.datasets, file = "/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/mrgn_v_gmac_v_mrpc_10k_datasets_c_kc_all_mods.RData")

#perform GMAC inferences
#GMAC
print("Running GMAC.....")
num.trios=length(small.datasets)
small.datasets.red=lapply(small.datasets, function(x)x[,1:3])
small.datasets.mat=do.call("cbind", small.datasets.red)
snp.idx=seq(1, dim(small.datasets.mat)[2]-2, 3)
snp.dat.cis=small.datasets.mat[,snp.idx]
exp.dat=small.datasets.mat[,-snp.idx]
trios.idx=cbind(c(1:num.trios), matrix(c(1:dim(exp.dat)[2]), nrow = num.trios, ncol = 2, byrow = T))


dim(t(exp.dat))
dim(t(snp.dat.cis))
dim(t(pc.matrix))

input.list=list(kc = t(kc), cov.pool=t(pc.matrix), snp.dat.cis=t(snp.dat.cis), exp.dat=t(exp.dat), trio.indexes=trios.idx,
                nperm=10000, use.nominal.p=TRUE)

cl <- makeCluster(10)

output.cis.med <- gmac(cl=input.list$clusters, known.conf = input.list$kc, cov.pool = input.list$cov.pool,
                       exp.dat = input.list$exp.dat, snp.dat.cis = input.list$snp.dat.cis,
                       trios.idx = input.list$trio.indexes, nperm = input.list$nperm,
                       nominal.p = input.list$use.nominal.p)

output.trans.med <- gmac(cl=input.list$clusters, known.conf = input.list$kc, cov.pool = input.list$cov.pool,
                       exp.dat = input.list$exp.dat, snp.dat.cis = input.list$snp.dat.cis,
                       trios.idx = input.list$trio.indexes[,c(1,3,2)], nperm = input.list$nperm,
                       nominal.p = input.list$use.nominal.p)
stopCluster(cl)

#reorganize output for saving
#cis results
out.table.cis=cbind.data.frame(output.cis.med[[1]], output.cis.med[[2]])
colnames(out.table.cis)=c(paste0('pval_', colnames(output.cis.med[[1]])),
                      paste0('effect_change_', colnames(output.cis.med[[2]])))

out.list.cis=list(out.table, input.list, output.cis.med[[3]])
names(out.list.cis)=c("output.table", "input.list", "cov.indicator.list")

#trans results
out.table.trans=cbind.data.frame(output.trans.med[[1]], output.trans.med[[2]])
colnames(out.table.cis)=c(paste0('pval_', colnames(output.trans.med[[1]])),
                          paste0('effect_change_', colnames(output.trans.med[[2]])))

out.list.trans=list(out.table.trans, input.list, output.trans.med[[3]])
names(out.list.trans)=c("output.table", "input.list", "cov.indicator.list")



print("Done!...Saving...")

save(out.list.trans, file = "/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/gmac_10k_trans_results_c_kc_all_mods.RData")
save(out.list.cis, file = "/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/gmac_10k_cis_results_c_kc_all_mods.RData")


#preform regressions and classify model types
#MRGN
print("Running MRGN")
reg.res=NULL
inf.mods=NULL
#regression
reg.res=sapply(small.datasets, infer.trio)
#model class
inf.mods=apply(reg.res,2, class.vec)

#get estimate of time to compute each trio
# for(i in 1:length(small.datasets)){
#   start.time=Sys.time()
#   z=infer.trio(small.datasets[[i]])
#   x=class.vec(z)
#   end.time=Sys.time()
#   params$Time.to.compute.mrgn[i]=difftime(end.time, start.time, units = 'mins')
# }
print("Done!...Saving results...")

#save
save(reg.res, file = "/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/mrgn_10k_regres_results_c_kc_all_mods.RData")
save(inf.mods, file = "/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/mrgn_10k_inf_results_c_kc_all_mods.RData")


#MRPC
print("Running MRPC...")
#----begin analysis-----#
#truth for M0
#V1-->T1
Truth.M0 <- MRPCtruth$M0
Adj.M0<- as(Truth.M0,"matrix")
#V1-->T2
Adj.M01 <- matrix(0,nrow=3,ncol = 3)
rownames(Adj.M01) <- colnames(Adj.M01) <- colnames(Adj.M0)
Adj.M01[1,3] <- 1

#truth for M1
#V1-->T1-->T2
Truth.M1 <- MRPCtruth$M1
Adj.M1<- as(Truth.M1,"matrix")
#V1-->T2-->T1
Adj.M11 <- matrix(0,nrow=3,ncol = 3)
rownames(Adj.M11) <- colnames(Adj.M11) <- colnames(Adj.M1)
Adj.M11[1,3] <- 1
Adj.M11[3,2] <- 1

#truth for M2
#V1-->T1<--T2
Truth.M2 <- MRPCtruth$M2
Adj.M2<- as(Truth.M2,"matrix")
#V1-->T2<--T1
Adj.M21 <- matrix(0,nrow=3,ncol = 3)
rownames(Adj.M21) <-colnames(Adj.M21) <-colnames(Adj.M2)
Adj.M21[1,3] <- 1
Adj.M21[2,3] <- 1
#truth for M3
#V1-->T1, V1-->T2
Truth.M3 <- MRPCtruth$M3
Adj.M3 <- as(Truth.M3,"matrix")

#truth for M4
#V1-->T1, V1-->T2, T1--T2
Truth.M4 <- MRPCtruth$M4
Adj.M4 <- as(Truth.M4,"matrix")

############################################3
#a function to apply mrpc in parallel
apply.mrpc=function(X){


  n <- nrow (X)
  V <- colnames(X)     # Column names
  print(n)
  # Classical correlation
  suffStat <- list(C = cor(X, use = "complete.obs"),
                   n = n)

  MRPC.fit.FDR <- MRPC(X,
                       suffStat,
                       GV = 1,
                       FDR = 0.05,
                       alpha=0.01,
                       indepTest = 'gaussCItest',
                       labels = V,
                       FDRcontrol = "ADDIS",
                       verbose = FALSE)


  #plot(MRPC.fit_FDR)
  Adj.infe1 <- as( MRPC.fit.FDR@graph,"matrix")
  print(Adj.infe1)
  Adj.infe <- Adj.infe1[1:3,1:3] #only consider snp, cis, trans
  colnames(Adj.infe) <- rownames(Adj.infe) <- colnames(Adj.M01)

  #use adj to determine model
  model="Other"
  if(identical(Adj.M0,Adj.infe)){model="M0.1"}
  if(identical(Adj.M01,Adj.infe)){model="M0.2"}
  if(identical(Adj.M1,Adj.infe)){model="M1.1"}
  if(identical(Adj.M11,Adj.infe)){model="M1.2"}
  if(identical(Adj.M2,Adj.infe)){model="M2.1"}
  if(identical(Adj.M21,Adj.infe)){model="M2.2"}
  if(identical(Adj.M3,Adj.infe)){model="M3"}
  if(identical(Adj.M4,Adj.infe)){model="M4"}
  #handle others

  return(list(model=model, Adj=Adj.infe))

}

# #obtain a subset of the confounders used to generate the model
# #represents the case we dont have the exact correct model
# #compare MRPC and MRGN
# smsmall.datasets=lapply(small.datasets, function(x){
#   if((dim(x)[2]-3)>15){
#     y=x[,c(1:3, sample(4:dim(x)[2], round(runif(1, 1, 15))))]
#     return(y)
#   }else{
#     return(x)
#   }
# })

#save(smsmall.datasets, file = "/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/mrpc.small.datasets.1k.RData")


#run the MRPC inference
#mrpc.infer=lapply(small.datasets, apply.mrpc)

#get estimate of time to infer each trio
mrpc.infer.list=list()
for(i in 1:length(small.datasets)){
  start.time=Sys.time()
  mrpc.infer.list[[i]]=apply.mrpc(small.datasets[[i]])
  end.time=Sys.time()
  params$Time.to.compute.mrpc[i]=difftime(end.time, start.time, units = 'mins')
}


print("Done!...Saving results...")

#save
save(mrpc.infer.list, file = "/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/mrpc_10k_small_results_c_kc_all_mods.RData")
save(params, file = "/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/mrpc_v_mrgn_v_gmac_10k_params_c_kc_all_mods.RData")

























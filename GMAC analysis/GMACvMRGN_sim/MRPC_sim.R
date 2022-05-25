
#Simulate All Data and Parameters:
library("GMAC", lib="/mnt/ceph/jarredk/Rpackages")
#cl <- makeCluster(10)
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
library("MRGN",lib="/mnt/ceph/jarredk/Rpackages")
library("qvalue", lib="/mnt/ceph/jarredk/Rpackages")

small.datasets=loadRData(file = "/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/simulated_data/mrgn_v_gmac_v_mrpc_10k_datasets_c_kc_all_mods.RData")
kc = t(as.matrix(output.WB$input.list$known.conf))

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
  #params$Time.to.compute.mrpc[i]=difftime(end.time, start.time, units = 'mins')
}


print("Done!...Saving results...")

#save
save(mrpc.infer.list, file = "/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/simulated_data/mrpc_10k_small_results_c_kc_all_mods.RData")
#save(params, file = "/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/simulated_data/mrpc_v_mrgn_v_gmac_10k_params_c_kc_all_mods.RData")



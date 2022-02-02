
# a script to simulate multiple types of trios and their indicator vector under the 
#MRPC updated framework to determine if the current set of cases is enough to identify
#each of the model types M0,M1,M2,M3,M4

#load necessary packages
#relevant packages
# library(MRPC)
# library(psych)
# library(prodlim)
library('MRPC', lib="/mnt/ceph/jarredk/Rpackages")
library('qvalue', lib="/mnt/ceph/jarredk/Rpackages")
library('psych', lib="/mnt/ceph/jarredk/Rpackages")
library('prodlim', lib="/mnt/ceph/jarredk/Rpackages")

# source("/mnt/ceph/jarredk/MRPC_UPDATE/MRPCgeneral.R")
# source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
source("/mnt/ceph/jarredk/MRPC_UPDATE/MRPC_update_simulation_functions.R")


#set initial parameters for simulation
#sample sizes
n=c(50, 100, 500, 1000)
#noise in the data/SD of errors
noise=c(0.2, 0.5, 0.8, 1.2)
#frequency of the minor allele
minor.allele=c(0.02, 0.1, 0.2, 0.3)
#signal strength of edge
b1.1=c(0.2, 0.4, 0.6, 0.8)
b1.2=c(0.2, 0.4, 0.6, 0.8)

#expand to all get all combos of simulation conditions
model.params=expand.grid(n, noise, minor.allele, b1.1, b1.2)
colnames(model.params)=c("sample.size","SD", "minor.freq", "b1.1","b1.2")
#the number of simulations for each set of conditions
sims=100




####################################################################

##----------M0-----------

all.data=list()
#data.sets=list()
accuracy=NULL

#looping for all combos of simulation conditions
for(i in 1:dim(model.params)[1]){
  #all.sim=list()
  samples100=list()
  #run simulations
  for(j in 1:sims){
    #simulate trio under parameters
    samples100[[j]]=SimulateData(N=model.params$sample.size[i], 
                                 p=model.params$minor.freq[i], 
                                 model="model0", 
                                 b0.1=0, 
                                 b1.1=model.params$b1.1[i], 
                                 sd.1=model.params$SD[i])
    
  }
  
  all.data[[i]]=samples100
  
  print(paste0("finished sim ", i, " of ", dim(model.params)[1]))
  
}

reg.res=list()
inf.mods=list()

for(i in 1:length(all.data)){
  
  reg.res[[i]]=sapply(all.data[[i]], regress)
  inf.mods[[i]]=sapply(reg.res[[i]], class.vec)
  
}

#pre-allocate stats
mean.acc=NULL
mean.acc4=NULL
mean.acc3=NULL
mean.acc2=NULL
#average accruacy across sample size
for(i in 1:length(n)){mean.acc[i]=mean(accuracy[which(model.params$sample.size==n[i])])}

#average accruacy across error SD
for(i in 1:length(noise)){mean.acc2[i]=mean(accuracy[which(model.params$SD==noise[i])])}

#average accruacy across b1.1 coeff
for(i in 1:length(b1.1)){mean.acc3[i]=mean(accuracy[which(model.params$b1.1==b1.1[i])])}

#average accruacy across minor allele freq.
for(i in 1:length(minor.allele)){mean.acc4[i]=mean(accuracy[which(model.params$minor.freq==minor.allele[i])])}


par(mfrow=c(2,2))

plot(minor.allele, mean.acc4, type="b", pch=21, bg="black",
     main="Avg. Accruacy Across Minor Allele Freq")

plot(b1.1, mean.acc3, type="b", pch=21, bg="black",
     main="Avg. Accruacy Across Signal")

plot(noise, mean.acc2, type="b", pch=21, bg="black",
     main="Avg. Accruacy Across Error SD")

plot(n, mean.acc, type="b", pch=21, bg="black",
     main="Avg. Accruacy Across Sample Size")
par(mfrow=c(1,1))

#summary(as.factor(store.vectors[[1]]$inf.class))

##----------M1-----------





##----------M2-----------



##----------M3-----------





##----------M4-----------
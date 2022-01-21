
# a script to simulate multiple types of trios and their indicator vector under the 
#MRPC updated framework to determine if the current set of cases is enough to identify
#each of the model types M0,M1,M2,M3,M4

#load necessary packages

library('MRPC', lib="/mnt/ceph/jarredk/Rpackages")
library('qvalue', lib="/mnt/ceph/jarredk/Rpackages")

source("/mnt/ceph/jarredk/MRPC_UPDATE/MRPCgeneral.R")
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")


simulate

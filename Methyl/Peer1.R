

#A simple helper function for easier reading in of .Rdata files with a specific name
#================================================Helper_Function============================================

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#read in expression data
bsgs=read.csv(file = "/mnt/ceph/jarredk/Methyl/new_data_GE_R2_bsgs_delete.csv")
bsgs.Rdata=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/new_data_GE_R2_bsgs_delete.RData")
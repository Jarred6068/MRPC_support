
#run trio specific algorithm on the same 5 tissues used in the GMAC analysis 
source("/mnt/ceph/jarredk/MRPC_UPDATE/GTEx_tissues_analysis/triospecific_GTEx.R")
top5=c(1,6, 33, 40, 48)
tissues.vec=tissue.names[top5, 2]
#run and save all data
result=mrnet.analyze(GTEx.tissname=tissues.vec, save=TRUE)






























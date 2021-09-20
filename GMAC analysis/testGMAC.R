
# library('MRPC')
# library('GMAC')
path1="C:/Users/Bruin/Documents/GitHub/MRPC_support/GMAC analysis/simMod4wSex.R"
source(path1)

library(qvalue)

greps=10
num.trios=100
N=1000
p.sex=0.5

for(i in 1:greps){
  
  exp.dat=vector("list", length = num.trios)
  known.conf=c(sample(c(1, 2), size = N, replace = TRUE, prob = c(p.sex, 1-p.sex)))
  snp.dat=as.data.frame(matrix(0, nrow = N, ncol = num.trios))
  
  for(j in 1:num.trios){
    
    sim.data <- sim.mod4.conf(N = N, 
                              p = 0.45,
                              b0.1 = 0,
                              b1.1 = 1, 
                              b1.2 = 1,
                              b1.3 = 0.05, 
                              sd.1 = 2,
                              Sex.vec=known.conf)
    print(head(sim.data))
    
    exp.dat[[j]]=sim.data[,-c(1,2)]
    
    snp.dat[,j]=sim.data[, 1]
    
    
  }
  
  trio.idx=cbind.data.frame(c(1:num.trios), matrix(c(1:(2*num.trios)), nrow = num.trios, ncol = 2, byrow = T) )
  colnames(trio.idx)=c("V1","T1","T2")
  
  exp.dat=as.data.frame(exp.dat)
  #print(dim(exp.dat))
  #print(str(exp.dat))
  
  num.M4=run.MRPC(exp.dat, snp.dat, trio.idx)
  print("--------------ADDIS---------------")
  print(num.M4)
  
  output <- gmac(known.conf = t(as.matrix(known.conf)), cov.pool = NULL, 
                 exp.dat = t(exp.dat), snp.dat.cis = t(snp.dat), 
                 trios.idx = trio.idx, nperm = 50, nominal.p = TRUE)
  
  qvals=qvalue(output$pvals$Known_sel_pool, fdr.level = 0.1)
  
  num.inferred=sum(qvals$significant)
  print("---------------GMAC----------------")
  print(num.inferred)
  
  
}

















#read in postprocessing functions
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")






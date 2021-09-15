
library("MRPC")
library("GMAC")


#N = 10^3;p = 0.45;b0.1 = 0;b1.1 = 1;b1.2 = 1;b1.3 = 0.05;sd.1 = 1;p.sex = 0.5

sim.mod4.conf=function(N, p, b0.1, b1.1, b1.2, b1.3, sd.1, Sex.vec){
  
  Sex=Sex.vec
  V1 <- c(sample(c(0, 1, 2), size = N, replace = TRUE, prob = c((1 - p)^2, 2 * p * (1 - p), p^2)))
  
  T1.a <- SimulateData1P(N = N, P1 = V1, b0.1 = b0.1, b1.1 = b1.1, 
                         sd.1 = sd.1)
  T2.a <- SimulateData3P(N = N, P1 = V1, P2 = T1.a, P3 = Sex, b0.1 = b0.1, 
                         b1.1 = b1.1, b1.2 = b1.2, b1.3=b1.3, sd.1 = sd.1)
  T2.b <- SimulateData1P(N = N, P1 = V1, b0.1 = b0.1, b1.1 = b1.1, 
                         sd.1 = sd.1)
  T1.b <- SimulateData3P(N = N, P1 = V1, P2 = T2.b, P3 = Sex, b0.1 = b0.1, 
                         b1.1 = b1.1, b1.2 = b1.2, b1.3=b1.3, sd.1 = sd.1)
  coinToss <- rbinom(n = N, size = 1, prob = 0.5)
  T1 <- rep(0, N)
  T1[which(coinToss == 0)] <- T1.a[which(coinToss == 0)]
  T1[which(coinToss == 1)] <- T1.b[which(coinToss == 1)]
  T2 <- rep(0, N)
  T2[which(coinToss == 0)] <- T2.a[which(coinToss == 0)]
  T2[which(coinToss == 1)] <- T2.b[which(coinToss == 1)]
  
  return(data.frame(V1 = V1, Sex=Sex, T1 = T1, T2 = T2))
  
}




run.MRPC=function(exp.dat, snp.dat, trio.idx){
  
  
  #truth for M4
  #V1-->T1, V1-->T2, T1--T2
  Truth.M4 <- MRPCtruth$M4
  Adj.M4 <- as(Truth.M4,"matrix")
  
  M44.FDR <- c()
  
  
  for(i in 1:dim(trio.idx)[1]){
    
    trio=cbind.data.frame(snp.dat[,unlist(trio.idx[i,1])], exp.dat[,unlist(trio.idx[i,2:3])])
    colnames(trio)=c("SNP","Cis","Trans")
  
    n <- nrow (trio.idx)
    V <- colnames(trio.idx)     # Column names
  
    # Classical correlation
    suffStat <- list(C = cor(trio, use = "complete.obs"),
                   n = n)
    
    print(cor(trio, use = "complete.obs"))
    #run MRPC on TRIO 
    MRPC.fit.FDR.addis <- MRPC(trio,
                               suffStat,
                               GV = 1,
                               FDR = 0.05,
                               indepTest = 'gaussCItest',
                               labels = V,
                               FDRcontrol = "ADDIS",
                               verbose = FALSE)
  
    MRPC.fit.addis.table=as(MRPC.fit.FDR.addis@graph, "matrix")
    Adj.infe=MRPC.fit.addis.table[1:3,1:3]
    #print(Adj.infe)
    #print(Adj.M4)
    #print(identical(Adj.M4,Adj.infe))
    
    if(identical(Adj.M4,Adj.infe)){
    M44.FDR[i]<-i
    }
    
    
  }
  
  return(length(M44.FDR))
  
  
}



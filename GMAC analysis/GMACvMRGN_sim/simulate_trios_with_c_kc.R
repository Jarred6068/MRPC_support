





#################################################
simData1P=function (N, P1, b0.1, b1.1, sd.1, G=NULL, u=NULL, kc=NULL, delta.kc=NULL, gamma.u=NULL) {

  if(is.null(c(u,kc))){
    t <- rnorm(n = N, mean = b0.1 + b1.1 * P1, sd = sd.1)
  }else{
    if(is.null(kc)){
      U=G[,u]
      #gamma.u=runif(length(u), min = cs.range[1], max=cs.range[2])
      t <- rnorm(n = N, mean = b0.1 + b1.1 * P1 + U%*%gamma.u, sd = sd.1)
    }else{
      U=G[,u]
      #gamma.u=runif(length(u), min = cs.range[1], max=cs.range[2])
      t <- rnorm(n = N, mean = b0.1 + b1.1 * P1 + U%*%gamma.u + kc%*%delta.kc, sd = sd.1)
    }

  }

  return(t)
}



#################################################
simData2P=function (N, P1, P2, b0.1, b1.1, b1.2, sd.1, G=NULL, u=NULL, kc=NULL, delta.kc=NULL, gamma.u=NULL) {

  if(is.null(c(u,kc))){
    t <- rnorm(n = N, mean = b0.1 + b1.1 * P1 + b1.2 * P2, sd = sd.1)
  }else{
    if(is.null(kc)){
      U=G[,u]
      #gamma.u=runif(length(u), min = cs.range[1], max=cs.range[2])
      t <- rnorm(n = N, mean = b0.1 + b1.1 * P1 + b1.2 * P2 + U%*%gamma.u, sd = sd.1)
    }else{
      U=G[,u]
      #gamma.u=runif(length(u), min = cs.range[1], max=cs.range[2])
      t <- rnorm(n = N, mean = b0.1 + b1.1 * P1 + b1.2 * P2 + U%*%gamma.u + kc%*%delta.kc, sd = sd.1)
    }

  }

  return(t)
}



#################################################
# simData3P=function (N, P1, P2, P3, b0.1, b1.1, b1.2, b1.3, sd.1, cs.max) {
#
#   t <- rnorm(n = N, mean = b0.1 + b1.1 * P1 + b1.2 * P2 + b1.3 *
#                P3, sd = sd.1)
#   return(t)
# }
#################################################




simDataNP=function (N, b0.1, sd.1, G=NULL, u=NULL, kc=NULL, delta.kc=NULL, gamma.u = NULL) {


  if(is.null(c(u,kc))){
    t <- rnorm(n = N, mean = b0.1, sd = sd.1)
  }else{
    if(is.null(kc)){
      U=G[,u]
      #gamma.u=runif(length(u), min = cs.range[1], max=cs.range[2])
      t <- rnorm(n = N, mean = b0.1 + U%*%gamma.u, sd = sd.1)
    }else{
      U=G[,u]
      #gamma.u=runif(length(u), min = cs.range[1], max=cs.range[2])
      t <- rnorm(n = N, mean = b0.1 + U%*%gamma.u + kc%*%delta.kc, sd = sd.1)
    }

  }
  return(t)

}





#################################################

simData=function (theta, model, b0.1, b.snp, b.med, sd.1, G=NULL, u=NULL, kc=NULL){

  #generate confounders
  G.dim=dim(G)[2]
  #conf=simConf(n = N, means = conf.means, variances = conf.var)
  #uniformly sample confounders
  #upper.b=round(2*u)
  #sample.u=round(runif(1, 1, upper.b))
  #sample.b=round(runif(1, 2, upper.b))
  #sample.c=round(runif(1, 2, upper.b))
  if(length(u)==1){
    u.idx=sample(1:G.dim, u, replace = FALSE)
  }else{
    u.idx=u
  }

  #b=sample(1:conf.dim, sample.b)
  #c=sample(1:dim(conf)[2], sample.c)
  N=dim(G)[1]

  #confounder effects - allow 50% to be negative and 50% to be positive
  #different coefficients for each node U is a parent of
  ct1 = rbinom(u, 1, 0.5)
  ct2 = rbinom(u, 1, 0.5)
  ct3 = rbinom(dim(kc)[2],1,0.5)
  ct4 = rbinom(dim(kc)[2],1,0.5)
  w1 = runif(u, min = 0.15, max = 0.5)
  w1[which(ct1==1)]=-1*w1[which(ct1==1)]
  w2 = runif(u, min = 0.15, max = 0.5)
  w2[which(ct2==1)]=-1*w2[which(ct2==1)]
  w3 = runif(dim(kc)[2], min = 0.01, max = 0.1)
  w3[which(ct3==1)]=-1*w3[which(ct3==1)]
  w4 = runif(dim(kc)[2], min = 0.01, max = 0.1)
  w4[which(ct4==1)]=-1*w4[which(ct4==1)]

  V1 <- c(sample(c(0, 1, 2), size = N, replace = TRUE,
                 prob = c((1 -theta)^2, 2 * theta * (1 - theta), theta^2)))

  colnames(G)=paste0("U",c(1:G.dim))

  switch(model, model0 = {

    T1 <- simData1P(N = N, P1 = V1, b0.1 = b0.1, b1.1 = b.snp, sd.1 = sd.1, G = G,
                    u = u.idx, gamma.u = w1, kc = kc, delta.kc = w3)
    T2 <- simDataNP(N = N, b0.1 = b0.1, sd.1 = sd.1, G = G, u = u.idx,
                    gamma.u = w2, kc = kc, delta.kc = w4)

    U=cbind.data.frame(G[,u.idx])
    return(data.frame(V1 = V1, T1 = T1, T2 = T2, U, kc))

  }, model1 = {

    T1 <- simData1P(N = N, P1 = V1, b0.1 = b0.1, b1.1 = b.snp, sd.1 = sd.1,
                    G = G, u = u.idx, gamma.u = w1, kc = kc, delta.kc = w3)
    T2 <- simData1P(N = N, P1 = T1, b0.1 = b0.1, b1.1 = b.med, sd.1 = sd.1,
                    G = G, u = u.idx, gamma.u = w2, kc = kc, delta.kc = w4)
    U=cbind.data.frame(G[,u.idx])
    return(data.frame(V1 = V1, T1 = T1, T2 = T2, U, kc))

  }, model2 = {

    T2 <- simDataNP(N = N, b0.1 = b0.1, sd.1 = sd.1, G = G, u = u.idx,
                    gamma.u = w1, kc = kc, delta.kc = w3)
    T1 <- simData2P(N = N, P1 = V1, P2 = T2, b0.1 = b0.1, b1.1 = b.snp, b1.2 = b.med, sd.1 = sd.1,
                    G = G, u = u.idx, gamma.u=w2, kc = kc, delta.kc = w4)
    U=cbind.data.frame(G[,u.idx])
    return(data.frame(V1 = V1, T1 = T1, T2 = T2, U, kc))

  }, model3 = {

    T1 <- simData1P(N = N, P1 = V1, b0.1 = b0.1, b1.1 = b.snp, sd.1 = sd.1,
                    G = G, u = u.idx, gamma.u=w1, kc = kc, delta.kc = w3)
    T2 <- simData1P(N = N, P1 = V1, b0.1 = b0.1, b1.1 = b.snp, sd.1 = sd.1,
                    G = G, u = u.idx, gamma.u=w2, kc = kc, delta.kc = w4)

    U=cbind.data.frame(G[,u.idx])
    return(data.frame(V1 = V1, T1 = T1, T2 = T2, U, kc))

  }, model4 = {

    T1.a <- simData1P(N = N, P1 = V1, b0.1 = b0.1, b1.1 = b.snp, sd.1 = sd.1,
                      G = G, u = u.idx, gamma.u=w1, kc = kc, delta.kc = w3)
    T2.a <- simData2P(N = N, P1 = V1, P2 = T1.a, b0.1 = b0.1, b1.1 = b.snp, b1.2 = b.med, sd.1 = sd.1,
                      G = G, u = u.idx, gamma.u=w2, kc = kc, delta.kc = w4)
    T2.b <- simData1P(N = N, P1 = V1, b0.1 = b0.1, b1.1 = b.snp, sd.1 = sd.1,
                      G = G, u = u.idx, gamma.u=w2, kc = kc, delta.kc = w4)
    T1.b <- simData2P(N = N, P1 = V1, P2 = T2.b, b0.1 = b0.1, b1.1 = b.snp, b1.2 = b.med, sd.1 = sd.1,
                      G = G, u = u.idx, gamma.u=w1, kc = kc, delta.kc = w3)

    coinToss <- rbinom(n = N, size = 1, prob = 0.5)
    T1 <- rep(0, N)
    T1[which(coinToss == 0)] <- T1.a[which(coinToss == 0)]
    T1[which(coinToss == 1)] <- T1.b[which(coinToss == 1)]
    T2 <- rep(0, N)
    T2[which(coinToss == 0)] <- T2.a[which(coinToss == 0)]
    T2[which(coinToss == 1)] <- T2.b[which(coinToss == 1)]
    U=cbind.data.frame(G[,u.idx])
    return(data.frame(V1 = V1, T1 = T1, T2 = T2, U, kc))

  }, stop("Model not included or missing"))
}








#load in normalized methylation data
loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

Mdata=loadRData( fileName = "/mnt/ceph/jarredk/Methyl/MethylData.RegressResids.Rdata")





calc.dens=function(data=NULL, mvec=NULL, sigvec=NULL, pvec=NULL, mix=NULL){
  
  dens.vec=NULL
  
  for(j in 1:mix){ dens.vec[j]=dnorm(data, mvec[j], sigvec[j])*pvec[j]}
  #print(dens.vec)
  return(sum(dens.vec))
}




EM=function(data=NULL, start.props=c(0.4, 0.6), start.mu=c(1, 19), start.sigmas=c(1,1), mixtures=2, conv.thresh=0.01){
  
  #initialize
  n=length(data)
  prop.vec=start.props
  mu.vec=start.mu
  sigma.vec=start.sigmas
  T.mat=as.data.frame(matrix(0, nrow = n, ncol = mixtures))
  
  
  prop.estimates=list()
  #prop.estimates[[1]]=start.props
  means.estimates=list()
  #means.estimates[[1]]=start.mu
  sigma.estimates=list()
  #sigma.estimates[[1]]=sigma.vec
  
  k=1
  t=1
  while(abs(t)>conv.thresh){
    
    
    #E-step
    for(i in 1:n){
      
      for(j in 1:mixtures){
        
        T.ij=(dnorm(data[i], mu.vec[j], sigma.vec[j])*prop.vec[j])/
          calc.dens(data=data[i], mvec=mu.vec, sigvec=sigma.vec, pvec=prop.vec, mix=mixtures)
        
        T.mat[i,j]=T.ij
          
      }
      
    }
    
    #print(T.mat)
    
    #M-step
    
    prop.estimates[[k]]=colSums(T.mat)/n
    means.estimates[[k]]=colSums(T.mat*data)/colSums(T.mat)
    sigma.estimates[[k]]=sqrt(colSums(T.mat*(data-means.estimates[[k]])^2)/colSums(T.mat))
    
    t=prop.estimates[[k]][1]-prop.vec[1]
    prop.vec=prop.estimates[[k]]
    mu.vec=means.estimates[[k]]
    sigma.vec=sigma.estimates[[k]]
    
    print(rbind(prop.vec, mu.vec, sigma.vec))
    #print(mu.vec)
    #print(sigma.vec)
    
    #print(t)
    k=k+1
    
  }
  
  #print(prop.estimates)
  return(T.mat)

}


#test algorithm

test.data=c(rnorm(2000, 0, 1), rnorm(8000, 20, 1))
hist(test.data)

T.mat1=EM(data = test.data)



#create fake data
a=c(rnorm(1000, 0.1, 1.5))
b=c(rnorm(200, 0, 1), rnorm(800, 20, 3))
c=c(rnorm(200, 0, 1), rnorm(200, 10, 1), rnorm(600, 20, 0.2))
fake1=cbind.data.frame(a,b,c)
colnames(fake1)=c("1G", "2G", "3G")


Mdata=fake1

#using Mclust

#install.packages("mclust", lib = "/mnt/ceph/jarredk/Rpackages", dependencies = T)

library('mclust', lib="/mnt/ceph/jarredk/Rpackages")


its=3
cp1=as.data.frame(matrix(0, nrow = its, ncol = 2))
cp2=as.data.frame(matrix(0, nrow = its, ncol = 2))
cp3=as.data.frame(matrix(0, nrow = its, ncol = 2))

colnames(cp1)=c("Loglik E", "Loglik V")
colnames(cp2)=c("Loglik E", "Loglik V")
colnames(cp3)=c("Loglik E", "Loglik V")

W.GLRT=as.data.frame(matrix(0, nrow = its, ncol = 4))
colnames(W.GLRT)=c("1v2 E", "1v2 V", "1v3 E", "1v3 V")


for(i in 1:its){
  
  hist(Mdata[,i], main = paste("hist of ", colnames(Mdata[,i])) )
  
  cp1[i,] = mclustLoglik(mclustBIC(Mdata[,i], G=1))
  cp2[i,] = mclustLoglik(mclustBIC(Mdata[,i], G=2))
  cp3[i,] = mclustLoglik(mclustBIC(Mdata[,i], G=3))
  
  W.GLRT[i,1:2]=-2*(cp1[i,]-cp2[i,])
  W.GLRT[i,3:4]=-2*(cp1[i,]-cp3[i,])
  
}

P.values=as.data.frame(matrix(0, nrow = its, ncol = 4))
colnames(P.values)=c("1v2 E", "1v2 V", "1v3 E", "1v3 V")

dfs=c(2,4,4,7)

for( i in 1:dim(W.GLRT)[2]){
  P.values[,i]=1-pchisq(W.GLRT[,i], dfs[i])
}







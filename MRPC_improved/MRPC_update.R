
#Data:

tissue="WholeBlood"

#read in badsha peerdata.v8 file
Peerdata.V8=loadRData(paste0("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/PEER_Files_V8/Peerdata.",tissue,".V8.RData"))

#preform PCA in same method as badsha and retain eigen values
Peerdata2 <- t(Peerdata.V8)
Peerdata3 <- Peerdata2[,apply(Peerdata2, 2, var, na.rm=TRUE) != 0]
# pca matrix
PCs<- prcomp(Peerdata3, scale=TRUE)$x

file1=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",
            tissue,"_AllPC/data.snp.cis.trans.final.",
            tissue,".V8.unique.snps.RData", sep = "")

trios=loadRData(fileName=file1)


#MRPC simulated data for model 1

library(MRPC)
data <- simu_data_M1 # load data for model 1
n <- nrow(data)      # Number of row
V <- colnames(data)  # Column names


library(qvalue, lib="/mnt/ceph/jarredk/Rpackages")
library(onlineFDR, lib="/mnt/ceph/jarredk/Rpackages")

ModelSkel=function(df=NULL, 
                   GV=NULL,
                   conf.pool=NULL, 
                   filter_fdr=0.05, 
                   labels=NULL,
                   FDRcontrol=c("LOND","ADDIS"), 
                   indepTest = c("gaussCItest", "disCItest", "citest"),
                   alpha=0.05)
  {
  
  p=dim(df)[2]
  
  cormat=cor(df, use="complete.obs")
  
  skel.adj=matrix(rep(0, p*p), nrow = p, ncol = p)
  colnames(skel.adj)=labels
  row.names(skel.adj)=labels
  
  cor.p=corr.test(data, data, use = "pairwise.complete.obs", adjust = "fdr", alpha = alpha)$p
  cor.p=cor.p[upper.tri(cor.p)]
  
  skel.adj[upper.tri(skel.adj)]=cor.p<alpha
  
  for(i in 1:GV){
    for(j in 1:(p-GV)){
      
      if(GV<=1){data.x=data}else{data.x=data[,-which(c(1:GV)!=i)]}
      
      out=as.data.frame(summary(lm(data[,j+GV]~., data=data.x[,-(j+GV)]))$coefficients)$`Pr(>|t|)`[-1]
      print(out<alpha)
    
      skel.adj[i,j+GV]=ifelse(out[1]<alpha, TRUE, FALSE)
    
    }
  }
  
  
  print(skel.adj)
}

ModelSkel(df=data, GV=1, labels=colnames(data))
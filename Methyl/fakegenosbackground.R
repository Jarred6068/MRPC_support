

#construct a fake matrix of 4,000,000 SNPs/genotypes for artificial trio construction

fake_genos=function(csize=NULL, org_genos, expmat=NULL){
  # percent.type=as.data.frame(matrix(0, nrow = dim(org_genos)[2], ncol = 4))
  # 
  # for(i in 1:dim(org_genos)[2]){ 
  #   
  #   org_genos[,i]=as.factor(org_genos[,i]) 
  #   percent.type[i,]=summary(org_genos[,i])/sum(summary(org_genos[,i]))
  #   
  # }
  
  #avg.per=colMeans(percent.type)
  
  genos=sample(c(0,1,2), csize*dim(expmat)[1], replace = TRUE)
  
  genos.mat=matrix(genos, nrow = dim(expmat)[1], ncol = csize, byrow = T)
  
  return(genos.mat)
  
}

genos.matBIG=fake_genos(4000000, ex.genotypes, express.aligned)
save(genos.matBIG, file = "/mnt/ceph/jarredk/Methyl/fakegenosBIG.Rdata")

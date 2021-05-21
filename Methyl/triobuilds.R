

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#Example assembly of trios script

#extract the two residual matrices and gene names:
ex.genotypes=read.table(file="/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/BrainCerebellarHemisphere.snp.data.V8.txt")
methyl.resids=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/MethylData.RegressResids.Rdata")
express.resids=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/ExpressData/Peerdata.bsgs.RData")
ex.genotypes=as.data.frame(ex.genotypes[-1,])

rn1=row.names(express.resids)[-1]
rn2=row.names(methyl.resids)[-1]

express.genenames=as.character(express.resids[1,])
express.resids2=as.data.frame(apply(express.resids[-1,],2,as.numeric))
row.names(express.resids2)=rn1

#extract just the bsgs_# part from row names in methyl matrix for alignment with expression
rn2.new=NULL
for(i in 1:length(rn2)){rn2.new[i]=strsplit(rn2[i], "/")[[1]][1]}

methyl.genenames=as.character(methyl.resids[1,])
methyl.resids2=as.data.frame(apply(methyl.resids[-1,],2,as.numeric))
row.names(methyl.resids2)=rn2.new

#alignment of subjects in both matrices

matchrows=match(rn2.new, rn1) 
express.aligned=express.resids2[matchrows,]
dim(express.aligned)
dim(methyl.resids2)

#parse the methylation genenames into a list of names vectors
methyl.gn.parsed=strsplit(methyl.genenames, ";")

#a function to match gene names from expression to the parsed list of methylation probe genes
#to pre-index the trio construction

qc=function(X, V){ X==V }

match.gn.parsed=function(egn, mgnp){
  
  indexlist1=vector("list", length(egn))
  vec1=c(1:length(mgnp))
  vec11=NULL
  
  for(i in 1:length(egn)){
    
    loglist1=lapply(mgnp, qc, egn[i])
      
    vec11=vec1[unlist(lapply(loglist1, any))]
    
    indexlist1[[i]]=vec11
    
  }
  
  return(indexlist1)
  
}



#construct a fake matrix of genotypes for artifical trio construction

triobuildlist=match.gn.parsed(express.genenames, methyl.gn.parsed)
triobuildlist[[1:5]]

save(triobuildlist, file = "/mnt/ceph/jarredk/Methyl/triobuildlist.Rdata")

percent.type=as.data.frame(matrix(0, nrow = dim(ex.genotypes)[2], ncol = 4))

for(i in 1:dim(ex.genotypes)[2]){ 
  
  ex.genotypes[,i]=as.factor(ex.genotypes[,i]) 
  percent.type[i,]=summary(ex.genotypes[,i])/sum(summary(ex.genotypes[,i]))
  
}

avg.per=colMeans(percent.type)
genos.mat=as.data.frame(matrix(0, nrow = dim(express.aligned)[1], ncol = length(unique(express.genenames)) ))

for(i in 1:dim(genos.mat)[2]){
  
  genos.mat[,i]=sample(c(0,1,2,NA), dim(genos.mat)[1], replace = TRUE, prob = avg.per)

}



#construct all possible trios:

























































































































































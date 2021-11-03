

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#Example assembly of trios script

#extract the two residual matrices and gene names:
#ex.genotypes=read.table(file="/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/BrainCerebellarHemisphere.snp.data.V8.txt")
methyl.resids=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/MethylData.RegressResids.Rdata")
express.resids=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/ExpressData/Peerdata.bsgs.RData")
ex.genotypes=as.data.frame(ex.genotypes[-1,])

rn1=row.names(express.resids)[-1]
rn2=row.names(methyl.resids)

express.genenames=as.character(express.resids[1,])
express.resids2=as.data.frame(apply(express.resids[-1,],2,as.numeric))
row.names(express.resids2)=rn1

#extract just the bsgs_# part from row names in methyl matrix for alignment with expression
rn2.new=NULL
for(i in 1:length(rn2)){rn2.new[i]=strsplit(rn2[i], "/")[[1]][1]}

methyl.genenames=loadRData(fileName ="/mnt/ceph/jarredk/Methyl/MethylData/MethylData.RegressResids.UCSC_RefGene_Name.Rdata")
#methyl.resids2=as.data.frame(apply(methyl.resids[-1,],2,as.numeric))
row.names(methyl.resids2)=rn2.new

#alignment of subjects in both matrices

matchrows=match(rn2.new, rn1) 
express.aligned=express.resids2[matchrows,]
dim(express.aligned)
dim(methyl.resids2)

save(methyl.resids2, file = "/mnt/ceph/jarredk/Methyl/MethylData.methylresids2.Rdata")
save(express.aligned, file = "/mnt/ceph/jarredk/Methyl/express.resids.Rdata")

#parse the methylation genenames into a list of names vectors
methyl.gn.parsed=strsplit(methyl.genenames, ";")

#a function to match gene names from expression to the parsed list of methylation probe genes
#to pre-index the trio construction:

#this function is slow and can take some time (>10 hrs) to run

qc=function(X, V){ X==V }

match.gn.parsed=function(egn, mgnp){
  
  indexlist1=vector("list", length(egn))
  vec1=c(1:length(mgnp))
  vec11=NULL
  
  for(i in 1:length(egn)){
    
    loglist1=lapply(mgnp, qc, egn[i])
      
    vec11=vec1[unlist(lapply(loglist1, any))]
    
    indexlist1[[i]]=vec11
    print(indexlist1[[i]])
    
  }
  
  return(indexlist1)
  
}

#a list of the methylation probes matched to each gene in express.genenames
triobuildlist=match.gn.parsed(express.genenames, methyl.gn.parsed)



#save(triobuildlist, file = "/mnt/ceph/jarredk/Methyl/triobuildlist.Rdata")
#====================================================================================
#form trios
#express.aligned=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/triobuildlist.Rdata")
#methyl.resids2=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/triobuildlist.Rdata")
#triobuildlist=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/triobuildlist.Rdata")
load(file="/mnt/ceph/jarredk/Methyl/Wksp1.Rdata")

#construct a fake matrix of 4,000,000 SNPs/genotypes for artificial trio construction

# fake_genos=function(csize=NULL, org_genos, expmat=NULL, expnames=NULL){
#   percent.type=as.data.frame(matrix(0, nrow = dim(org_genos)[2], ncol = 4))
# 
#   for(i in 1:dim(org_genos)[2]){ 
#   
#     org_genos[,i]=as.factor(org_genos[,i]) 
#     percent.type[i,]=summary(org_genos[,i])/sum(summary(org_genos[,i]))
#   
#   }
# 
#   avg.per=colMeans(percent.type)
#   genos.mat=as.data.frame(matrix(0, nrow = csize, ncol = dim(expmat)[1] ))
# 
#   
#     genos=sample(c(0,1,2,NA), csize*dim(genos.mat)[2], replace = TRUE, prob = avg.per)
#     
#     genos.mat=matrix(genos, nrow = dim(genos.mat)[2], ncol = csize, byrow = T)
#     
#     return(genos.mat)
# 
# }

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

#genos.matBIG=fake_genos(4000000, ex.genotypes, express.aligned)
#save(genos.matBIG, file = "/mnt/ceph/jarredk/Methyl/fakegenosBIG.Rdata")

genos.mat=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/fakegenosBIG.Rdata")
genos.mat[1:10,1:5]

#library(gdata, lib="/mnt/ceph/jarredk/Rpackages")

#keep(express.genenames, express.aligned, express.resids2, 
#     methyl.genenames, methyl.gn.parsed, genos.mat, rn2.new, 
#     methyl.resids2, loadRData, triobuildlist, sure = TRUE)

#------------------plotting-to-check-normalization-----------------------

for(j in c(1,2,4,5,7:14)){
submat.methyl=methyl.resids2[, triobuildlist[[j]] ]

  for(i in 1:dim(submat.methyl)[2]){
  
    png(paste("/mnt/ceph/jarredk/Methyl/NormCheckPlots/","pairplots_",express.genenames[j],"_",colnames(express.aligned)[j],
              "_", colnames(submat.methyl)[i], ".png", sep = ""))
    plot(express.aligned[,j], submat.methyl[,i], 
        xlab = "Expression Residuals",
        ylab = "Methylation Residuals",
        main = paste("plot:", express.genenames[j],":",colnames(express.aligned)[j],
                     ":", colnames(submat.methyl)[i], sep = " "))
    abline(a=0, b=1, col="red")
  
    dev.off()
  
  }
}

for(j in c(1,2,4,5,7:14)){
  submat.methyl=methyl.resids2[, triobuildlist[[j]] ]
  
  for(i in 1:dim(submat.methyl)[2]){
    
    png(paste("/mnt/ceph/jarredk/Methyl/ClusCheckPlots/","plots", colnames(submat.methyl)[i], ".png", sep = ""))
    hist(submat.methyl[,i], 
         breaks = 20,
         xlab = "Post-normalization Methylation Residuals",
         main = paste("Histogram:",colnames(submat.methyl)[i], sep = " "))
    
    dev.off()
    
  }
}

#------------------------------------------------------------------------

#construct all trios:

start.time=Sys.time()

build.trios=function(emat=NULL, mmat=NULL, tbl=NULL, enames=NULL){
  #emat -- the expression data matrix
  #mmat -- the methylation data matrix
  #tbl -- the triobuildlist (output from match.gn.parsed)
  #enames -- the gene names that correspond to the columns of the 
  #           expression matrix
  
  trio.full.list=vector("list", length = length(enames))
  names(trio.full.list)=enames
  checker=integer(0)
  
  for(i in 1:length(enames)){
    print(i)
    
    if(identical(checker, tbl[[i]])){
      
      trio.full.list[[i]]=NA
      
    }else{
      
      #print(mmat[,tbl[[i]] ])
      submatM=as.matrix(mmat[, tbl[[i]] ])
      trio.gl=vector("list", length = length(tbl[[i]]) )
      
      for(j in 1:length(tbl[[i]])){
        combmat=cbind.data.frame(emat[,i], submatM[, j])
        colnames(combmat)=c(paste0(enames[i],";",colnames(emat)[i]), colnames(submatM)[j])
        row.names(combmat)=row.names(submatM)
        #print(combmat)
        trio.gl[[j]]=combmat
      }
      
      trio.full.list[[i]]=trio.gl
      
    }
    
    
  }
  
  return(trio.full.list)
  
}

trios=build.trios(emat = express.aligned, mmat = methyl.resids2, tbl = triobuildlist, enames = express.genenames)

end.time=Sys.time()

print(end.time-start.time)

save(trios, file = "/mnt/ceph/jarredk/Methyl/completedTrios1.Rdata")



















































































































































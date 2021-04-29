

#A simple helper function for easier reading in of .Rdata files with a specific name
#================================================Helper_Function============================================

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#read in expression data
bsgs=read.csv(file = "/mnt/ceph/jarredk/Methyl/ExpressData/new_data_GE_R2_bsgs_delete.csv")
bsgs=as.data.frame(bsgs[,-1])
bsgs.Rdata=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/ExpressData/new_data_GE_R2_bsgs_delete.RData")
geo.bsgs=read.csv(file = "/mnt/ceph/jarredk/Methyl/ExpressData/BSGS_GEO_Accession_data2.csv")

#methyl=read.csv(file = "/mnt/ceph/jarredk/Methyl/new_data_M_average_delete.csv")

#remove p-values
pcols=seq(4, dim(bsgs)[2], 2)
#separate the pcols and gene express cols
p.bsgs=t(bsgs[,pcols])
colnames(p.bsgs)=bsgs[,1]

#use the selection described in brisbane paper

percent.sig=NULL
idx=c(1:dim(p.bsgs)[2])

for(i in 1:dim(p.bsgs)[2]){
  
  percent.sig[i]=sum(p.bsgs[,i]<0.05)/length(p.bsgs[,i])
  
}
cols.to.keep=idx[percent.sig>=0.1]

#extract gene expression values
bsgs=bsgs[,-pcols]
#transpose so subjects are rows
tbsgs=t(bsgs)
Ref_id=tbsgs[1,]
Genes=tbsgs[2,]
tbsgs=tbsgs[-c(1,2),]
#rename cols and rows
subjects=row.names(tbsgs)
#convert character to numeric
tbsgs=as.data.frame(apply(tbsgs, 2, as.numeric))
colnames(tbsgs)=Ref_id
row.names(tbsgs)=subjects

tbsgs=tbsgs[, cols.to.keep]

tbsgs[1:10,1:10]

#matching up covariates with rows of bsgs data
covariates=na.omit(geo.bsgs[,3:6])
cov.id=covariates[,4]
idx2=NULL

for(i in 1:length(cov.id)){
  
  idx2[i]=match(subjects[i], cov.id)
  
}

covariates.final=covariates[idx,]
covariates.final[1:10,]

#center and scale
library('preprocessCore',lib="/mnt/ceph/jarredk/Rpackages")

tbsgs.scaled.normal=normalize.quantiles(scale(tbsgs))

save(covariates.final, file = "/mnt/ceph/jarredk/Methyl/ExpressData/bsgs.cov.final.Rdata")
save(tbsgs, file = "/mnt/ceph/jarredk/Methyl/ExpressData/transposed.bsgsdata.Rdata")
save(tbsgs.scaled.normal, file = "/mnt/ceph/jarredk/Methyl/ExpressData/transposed.bsgs.scaled.normal.Rdata")

#histograms to check normality

for( i in 1:100){
  #pre-transformation:
  png(paste("/mnt/ceph/jarredk/Methyl/ExpressData/Histograms/","histogram_", colnames(tbsgs)[i],i, ".png", sep = ""))
  H1=hist(tbsgs[,i])
  H2=hist(tbsgs.scaled.normal[,i])
  par(mfrow=c(1,2))
  plot(H1)
  plot(H2)
  par(mfrow=c(1,1))
  
  dev.off()
  
  
}


for( i in 1:100){
  #pre-transformation:
  png(paste("/mnt/ceph/jarredk/Methyl/ExpressData/Pair_plots/","pairplots_", colnames(tbsgs)[i],i, ".png", sep = ""))
  plot(tbsgs.scaled.normal[,i], tbsgs[,i], 
       xlab = "Scaled and quantile normalized",
       ylab = "original value",
       main = paste(colnames(tbsgs)[i]))
  
  dev.off()
  
  
}

#========================================================PEER==============================================
#dataprep

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


#bsgs=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/new_data_GE_R2_bsgs_delete.RData")
tbsgs=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/ExpressData/transposed.bsgsdata.Rdata")
tbsgs.scaled.normal=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/ExpressData/transposed.bsgs.scaled.normal.Rdata")
covariates.final=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/ExpressData/bsgs.cov.final.Rdata")

#1
library(peer,lib="/mnt/ceph/jarredk/Rpackages")

#2 PEER model object
model = PEER()

#15 factors for N < 150
#30 factors for 150 ??? N < 250
#45 factors for 250 ??? N < 350
#60 factors for N ??? 350
if(dim(tbsgs)[1]<150){K=15}
if(dim(tbsgs)[1]>=150 && dim(tbsgs)[1]< 250){K=30}
if(dim(tbsgs)[1]>=250 && dim(tbsgs)[1]< 350){K=45}
if(dim(tbsgs)[1]>=350){K=60}


#infer K=10 hidden confounders or peer factor 
#3
print(K)
PEER_setNk(model,K)

#4 set the observed gene expression data
PEER_setPhenoMean(model, as.matrix(tbsgs.scaled.normal))
dim(PEER_getPhenoMean(model))
#5 covariates matrix
PEER_setCovariates(model, as.matrix(covariates.final))

# 6 perform the inference
PEER_update(model)

#residuals
residuals.bsgs= PEER_getResiduals(model)
dim(residuals.bsgs)
residuals.bsgs[1:5,1:5]

colnames(residuals.bsgs)<-colnames(tbsgs)
rownames(residuals.bsgs)<-rownames(tbsgs)
residuals.bsgs[1:5,1:5]
#save
Peerdata.bsgs=residuals.bsgs
Peerdata.bsgs[1:5,1:5]
save(Peerdata.bsgs, file="/mnt/ceph/jarredk/Methyl/ExpressData/Peerdata.bsgs.RData")
write.csv(Peerdata.bsgs, file="/mnt/ceph/jarredk/Methyl/ExpressData/Peerdata.bsgs.csv")

# png("/mnt/ceph/jarredk/Methyl/hist3.png")
# hist(tbsgs.log[,3])
# 
# dev.off()





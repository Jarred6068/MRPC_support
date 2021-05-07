

#A simple helper function for easier reading in of .Rdata files with a specific name
#================================================Helper_Function============================================

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#read in methylation data, bsgs/brisbane subject ID equivalents, and covariates
methyl_original=read.csv(file = "/mnt/ceph/jarredk/Methyl/new_data_M_average_delete.csv")
geo.bsgs=read.csv(file = "/mnt/ceph/jarredk/Methyl/ExpressData/BSGS_GEO_Accession_data2.csv")
geo_id_equivs=read.table(file = "/mnt/ceph/jarredk/Methyl/GSE53195_GSE56105_ID_equivalence.txt", header = T)
methyl_original[1:10,1:5]
head(geo.bsgs)
head(geo_id_equivs)

#save relevant columns, column names, and row names
UCSC_RefGene_Name=methyl_original$UCSC_RefGene_Name
ID_Ref=methyl_original$ID_REF
#remove meta information columns
methyl=methyl_original[,-c(1:3)]
methyl[1:10,1:5]
#partition data into p-vals and methylation sets
p.val.cols=seq(2, dim(methyl)[2], 2)
methyl.pvals=t(methyl[,p.val.cols])
methyl.pvals[1:10,1:5]
methyl.M=methyl[,-p.val.cols]
methyl.M[1:10,1:5]

#get and match the subject IDs 
subject_IDs=colnames(methyl.M)
idx=match(geo_id_equivs$GSE56105_ID, subject_IDs)
#keep subjects that match, and transpose so subjects are in rows
tmethyl.M=t(methyl.M[,idx])
tmethyl.M[1:10, 1:5]
#check
all.equal(row.names(tmethyl.M), geo_id_equivs$GSE56105_ID)
#removal of columns with more than 5 rows having detection p-value > 0.001
#or 11 or more individuals with na
pvals1=NULL
missing=NULL
col.idx=1:dim(methyl.pvals)[2]

for(i in 1:dim(methyl.pvals)[2]){
  
  pvec=methyl.pvals[,i]
  mvec=tmethyl.M[,i]
  pvals1[i]=sum(pvec>0.001)
  missing[i]=sum(is.na(mvec))
  
}
info1=cbind.data.frame(pvals1, missing, col.idx)
colnames(info1)=c("p","m","idx")


cols.to.remove=subset(info1, p>5 | m>11)$idx
print(cols.to.remove)


#remove
tmethyl.M=tmethyl.M[,-cols.to.remove]

#align covariates
#matching up covariates with subjects (i.e rows) of methyl data
covariates=na.omit(geo.bsgs[,3:6])
cov.id=covariates[,4]

  idx2=match(cov.id, geo_id_equivs$GSE53195_ID)
  
covariates.final=covariates[na.omit(idx2),]
covariates.final[1:10,]
  
tmethyl.M=tmethyl.M[na.omit(idx2),]
dim(tmethyl.M)

#rename columns and rows
colnames(tmethyl.M)=ID_Ref[-cols.to.remove]
row.names(tmethyl.M)=paste0(geo_id_equivs$GSE53195_ID[na.omit(idx2)], 
                            rep("/", length(geo_id_equivs$GSE56105_ID[na.omit(idx2)])), 
                            geo_id_equivs$GSE56105_ID[na.omit(idx2)])

tmethyl.M[1:10, 1:5]

save(tmethyl.M, file = "/mnt/ceph/jarredk/Methyl/Mdata.matched.final.Rdata")
save(covariates.final, file = "/mnt/ceph/jarredk/Methyl/Mdata.matched.covars.final.Rdata")

#==========================================================================================
#read in the data

tmethyl.M = loadRData(fileName = "/mnt/ceph/jarredk/Methyl/Mdata.matched.final.Rdata")
covariates.final = loadRData(fileName = "/mnt/ceph/jarredk/Methyl/Mdata.matched.covars.final.Rdata")

#probes are proportions i.e between 0 and 1
#preform logit transformation on each probe

LT=function(X){
  
  Z=log(X/(1-X))
  
  return(Z)
}


tmethyl.Mlt=apply(tmethyl.M, 2, LT)


for( i in 1:100){
  #pre-transformation:
  png(paste("/mnt/ceph/jarredk/Methyl/MHistograms/","histogram_", colnames(tmethyl.M)[i],i, ".png", sep = ""))
  H1=hist(tmethyl.M[,i], main = "Untransformed")
  H2=hist(tmethyl.Mlt[,i], main = "logit transformed")
  par(mfrow=c(1,2))
  plot(H1)
  plot(H2)
  par(mfrow=c(1,1))
  dev.off()
  
}


#Regress each logit transformed probe against the covariates as described by:
#McRae, Allan F., et al. "Contribution of genetic variation to transgenerational 
#inheritance of DNA methylation." Genome biology 15.5 (2014): 1-10.

#convert gender to a factor
covariates.final$Gender=as.factor(covariates.final$Gender)

#=========================================function to impute the response=======================================================

#runs interal to get.resids
impute.avg=function(data.mat=NULL){
  
  #locate missing data
  omitted=na.omit(data.mat$log.Y)
  locations=attr(omitted, "na.action")
  
  #extract regressors
  regressors=data.mat[locations, -1]
  
  group.avg=NULL
  #impute average of obs with regressor traits
  for(i in 1:length(locations)){
    
    x=subset(data.mat, Gender == regressors[i,1] & Age == regressors[i,2])[,1]
    print(x)
    group.avg[i]=mean(na.omit(x))
    
  }
  
  data.mat$log.Y[locations] = group.avg
  #return imputed obs vector
  return(data.mat$log.Y)
  
}


#=========================================function to get and plot the residuals:===============================================
get.resids=function(X, covars, impute.method="column_average", plot.pairs=TRUE, verbose=0){
  #--------------------------------------------------------------------------------------------
  #SYNTAX:
  #X -- the matrix of filtered and logit transformed methylation probes
  #covars -- matrix of covariates
  #impute.method -- one of "column_average", "FAMD", or "PCA"
  #plot.pairs -- logical indicating if first 100 probes should be plotted against the residuals
  #verbose -- integer in [0, 1]: setting to 1 outputs progress printing
  #--------------------------------------------------------------------------------------------
  
  
  library(missMDA)
  
  residuals.matrix=as.data.frame(matrix(0, nrow = dim(tmethyl.Mlt)[1], ncol = dim(tmethyl.Mlt)[2]))
  #preform a regression on each methylation probe against covariates:
  
  
  for(i in 1:dim(X)[2]){
    #bind column and covars into easy-to-use dataframe for lm
    data1=cbind.data.frame(X[,i], covars)
    colnames(data1)=c("log.Y", colnames(covars))
    
    #impute any missing values:
    if(sum(is.na(data1$log.Y))>0){
      print("found an NA", sum(is.na(data1$log.Y)))
      
      data1$log.Y=impute.avg(data1)
      
      print(sum(is.na(data1$log.Y)))
      
    }
    
    #fit model 
    print(head(data1))
    print(summary(data1))
    
    model=lm(log.Y ~ Gender*poly(Age, 2), data = data1)
    
    print(summary(model))
    #get residuals
    residuals.matrix[,i]=model$residuals
    
  }
  
  #plotting
  if(plot.pairs==TRUE){

    for( i in 1:100){
      #pre-transformation:
      png(paste("/mnt/ceph/jarredk/Methyl/ExpressData/Pair_plotsM/","pairplots_",colnames(X)[i],i, ".png", sep = ""))
      plot(X[,i], residuals.matrix[,i], 
           xlab = "Logit Transformed",
           ylab = "Regression Residuals",
           main = paste("plot:", colnames(X)[i], sep = " "))
      abline(0, 1, col="red")
      
      dev.off()
      
      
    }
  }
  
  colnames(residuals.matrix)=colnames(X)
  row.names(residuals.matrix)=row.names(X)
  #return the residual matrix (data adjusted for covariates)
  return(residuals.matrix)
    
}

#====================================================================================================================



































source("/mnt/ceph/jarredk/Reg_Net/AL_genetabV2.R")
#sort.gt2()
#sort.gt2(mediator="cis")
#sort.gt2(FDR="LOND",mediator="cis")
#sort.gt2(FDR="LOND",mediator="trans")

#create master df
df=master.create()
dim(df$lm1t1)
dim(df$lm1t2)
dim(df$am1t1)
dim(df$am1t2)

#counts and bins for 4 subsets
LTasM=binit(df$lm1t2, target=4)
LCasM=binit(df$lm1t1, target=3)

ATasM=binit(df$am1t2, target=4)
ACasM=binit(df$am1t1, target=3)


#Helper plotting function
plotit=function(data=NULL, name=NULL){
  
  png(paste("/mnt/ceph/jarredk/Reg_Net/Figures/rplot_", name  ,".png", sep = ""))
  
  H1=hist(data$counts[,2], breaks = 15)
  plot(H1, xlab = "# of Tissues Shared", ylab ="# of Genes",
       main = paste("Frequency of Shared Tissues for",name,sep = " "))
  
  dev.off()
  
}

#plot the histograms for 4 subsets
plotit(LTasM, "LTasM")
plotit(LCasM, "LCasM")
plotit(ATasM, "ATasM")
plotit(ACasM, "ACasM")


#creating some summary tables
library(knitr, lib="/mnt/ceph/jarredk/Rpackages")
#hist--tables
#===============================================================================
sort.bins=function(bins.table=NULL){
  
  G=sort(bins.table[,2], index.return=TRUE)
  bins.table=bins.table[G$ix,]
  
  return(bins.table)
}

#===============================================================================

kable(sort.bins(LTasM$bins), caption="LOND Trans Mediated")
kable(sort.bins(LCasM$bins), caption="LOND Cis Mediated")
kable(sort.bins(ATasM$bins), caption="ADDIS Trans Mediated")
kable(sort.bins(ACasM$bins), caption="ADDIS Cis Mediated")


#percentage of unqiue genes---tables
#ADDIS
AT.unq.genes=dim(ATasM$counts)[1]
AC.unq.genes=dim(ACasM$counts)[1]

total.unq.genes=AT.unq.genes+AC.unq.genes

percent.AT=AT.unq.genes/total.unq.genes
percent.AC=AC.unq.genes/total.unq.genes

#LOND
LT.unq.genes=dim(LTasM$counts)[1]
LC.unq.genes=dim(LCasM$counts)[1]

total.unq.genesL=LT.unq.genes+LC.unq.genes

percent.LT=LT.unq.genes/total.unq.genesL
percent.LC=LC.unq.genes/total.unq.genesL

#TABLES
T1=cbind.data.frame(c(AT.unq.genes, percent.AT), c(AC.unq.genes, percent.AC))
T2=cbind.data.frame(c(LT.unq.genes, percent.LT), c(LC.unq.genes, percent.LC))

colnames(T1)=c("Trans.Mediated","Cis.Mediated")
row.names(T1)=c("Total.Unique","% Unique")
colnames(T2)=c("Trans.Mediated","Cis.Mediated")
row.names(T2)=c("Total.Unique","% Unique")

kable(T1, caption="ADDIS Results")
kable(T2, caption="LOND Results")

#Helper function
match.rowN=function(df=NULL, target.vec=NULL){
  
  idx=rep(0, length(target.vec))
  
  for(i in 1:length(target.vec)){
    
    idx[i]=match(target.vec[i], df)
    
  }
  
  return(idx)
  
}

#percentage of cis and trans mediator gene types for ADDIS

IX=match.rowN(df=df$am1t1$Cis.Gene.ID, target.vec = ACasM$counts$Gene.stable.ID)
am1t1.unique=df$am1t1[IX,]

IX=match.rowN(df=df$am1t2$Trans.Gene.ID, target.vec = ATasM$counts$Gene.stable.ID)
am1t2.unique=df$am1t2[IX,]

am1t1.unique$cis.Gene.type=as.factor(am1t1.unique$cis.Gene.type)
am1t2.unique$trans.Gene.type=as.factor(am1t2.unique$trans.Gene.type)

percent.ACasM.types=summary(am1t1.unique$cis.Gene.type)/sum(summary(am1t1.unique$cis.Gene.type))
percent.ATasM.types=summary(am1t2.unique$trans.Gene.type)/sum(summary(am1t2.unique$trans.Gene.type))

kable(percent.ACasM.types, col.names="Percentage", caption="ADDIS Cis Mediator Gene Types")
kable(percent.ATasM.types, col.names="Percentage", caption="ADDIS Trans Mediator Gene Types")




# r=getOption("repos")
# install.packages(knitr, lib = "/mnt/ceph/jarredk/Rpackages", repos=r["CRAN"])



















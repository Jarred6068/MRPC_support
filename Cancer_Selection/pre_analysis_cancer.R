
#####################################################################

source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")

#read in TCGA Methylation data
TCGA.meth=loadRData("/mnt/ceph/jarredk/Cancer_Selection/TCGA.meth.RData")
TCGA.meta=read.table("/mnt/ceph/jarredk/Cancer_Selection/brca_tcga_clinical_data_clean_1.txt", sep="\t", header=T)

TCGA.gene.name=TCGA.meth$Gene_Symbol
TCGA.gene.coord=TCGA.meth$Genomic_Coordinate
TCGA.gene.chr=TCGA.meth$Chromosome

#extract the samples for which detailed cancer subtype is available:
samples1=colnames(TCGA.meth)[-which(is.na(TCGA.meta$Cancer.Type.Detailed))]
matched=match(samples1,colnames(TCGA.meth))
TCGA.meth.matched=TCGA.meth[,matched]
# TCGA.meth.matched.subtype=cbind.data.frame(TCGA.meth.matched[,1:4], 
#                                            Subtype=TCGA.meta$Cancer.Type.Detailed[-which(is.na(TCGA.meta$Cancer.Type.Detailed))],
#                                            TCGA.meth.matched[,-c(1:4)])

#read in BSGS Methylation data
BSGS.meth=loadRData("/mnt/ceph/jarredk/Methyl/MethylData/MethylData.RegressResids.Rdata")
BSGS.meta=read.csv("/mnt/ceph/jarredk/Methyl/GPL13534_HumanMethylation450_15017482_v.1.1.csv")


#Calculate the PCA matrix for TCGA methyl data
#look for potential clustering via cancer subtype

TCGA.pca=prcomp(t(TCGA.meth), retx=T)$x

png("/mnt/ceph/jarredk/Cancer_Selection/TCGA.pca.plot.png")

plot(TCGA.pca[,1], TCGA.pca[,2], color=as.factor(TCGA.meta$))


#####################################################################



#process stem cell data 
library('BiocManager',lib="/mnt/ceph/jarredk/Rpackages")
BiocManager::install("hgu133a.db", type = "source", lib="/mnt/ceph/jarredk/Rpackages")
BiocManager::install("hgu133acdf", lib="/mnt/ceph/jarredk/Rpackages")
BiocManager::install("AnnotationDbi", lib="/mnt/ceph/jarredk/Rpackages")
BiocManager::install("illuminaio", lib="/mnt/ceph/jarredk/Rpackages")



#####################################################################
#reading hESC expressiong data and normalizing via RMA
library('affy', lib="/mnt/ceph/jarredk/Rpackages")
library('GEOquery', lib="/mnt/ceph/jarredk/Rpackages")
setwd("/mnt/ceph/jarredk/Cancer_Selection/GSE7332_Stem_cell/")
#download files from GEO
getGEOSuppFiles("GSE7332")
setwd("/mnt/ceph/jarredk/Cancer_Selection/GSE7332_Stem_cell/GSE7332/")
#untar and unzip all .cel files
untar("GSE7332_RAW.tar", exdir = "data")
cels = list.files("data/", pattern = "CEL")
sapply(paste("data", cels, sep = "/"), gunzip)
cels = list.files("data/", pattern = "CEL")

#load platform info
library("org.Hs.eg.db", lib="/mnt/ceph/jarredk/Rpackages")
library("hgu133a.db", lib="/mnt/ceph/jarredk/Rpackages")
library("hgu133acdf", lib="/mnt/ceph/jarredk/Rpackages")
setwd("/mnt/ceph/jarredk/Cancer_Selection/GSE7332_Stem_cell/GSE7332/data/")

#preform RMA normalization
raw.data = ReadAffy(verbose = FALSE, filenames = cels, cdfname = "hgu133acdf")
data.rma.norm = rma(raw.data)
rma = exprs(data.rma.norm)
rma[1:5,1:5]
write.table(rma, file = "rma.txt", quote = FALSE, sep = "\t")

#combine rma with probe ids
tt = cbind.data.frame(row.names(rma), rma)
colnames(tt) = c("ProbeID", sub(".cel", "", colnames(rma), ignore.case = TRUE))
rownames(tt) = NULL
tt[1:5, 1:5]

#get annotations and gene names
library("AnnotationDbi",lib="/mnt/ceph/jarredk/Rpackages")
OUT <- select(hgu133a.db,keys= as.character(tt[,1]), columns=c("PROBEID","SYMBOL", "ENTREZID", "GENENAME"))

matched=match(tt$ProbeID, OUT$PROBEID)

#join gene info with expression info
tt.annot=cbind.data.frame(OUT[matched,],tt[,-1])

#save
write.table(tt.annot, file = "rma.annotated.txt", quote = FALSE, sep = "\t")

#####################################################################

#reading normal human breast tissue files


#total of 121 healthy human breast tissue samples:
#read in signal intensities matrix (methylated sites count, unmethylated sites ct, and pvalues)
setwd("/mnt/ceph/jarredk/Cancer_Selection/GSE101961_normal_breast/")
mat.sig=read.table("GSE101961_Matrix_signal.txt", header=T, sep="\t")

#parse colnames 
names.list=strsplit(colnames(mat.sig), split=".", fixed = TRUE)
#get the pval columns:
methyl.sig.pcols=lapply(names.list, FUN = function(x)  any(x=="Pval"))
mat.sig.pval.only=mat.sig[,c(1, which(methyl.sig.pcols==TRUE))]
dim(mat.sig.pval.only) #dimension should be probes x samples + ID column

#which columns are methylated sig cts -- allocate to separate matrix
methyl.sig.tf=lapply(names.list, FUN = function(x)  any(x=="Methylated"))
mat.sig.methyl=mat.sig[,c(1, which(methyl.sig.tf==TRUE))]
dim(mat.sig.methyl) #dimension should be probes x samples + ID column

#which columns are methylated sig cts -- same as above but for unmethylated
unmethyl.sig.tf=lapply(names.list, FUN = function(x)  any(x=="Unmethylated"))
mat.sig.unmethyl=mat.sig[,c(1,which(unmethyl.sig.tf==TRUE))]
dim(mat.sig.unmethyl)

#calculate % of methylation at each probe site for each sample

final.methyl=cbind.data.frame(ProbeID=mat.sig[,1], mat.sig.methyl[,-1]/(mat.sig.methyl[,-1]+mat.sig.unmethyl[,-1]))
#save
write.table(final.methyl, file="/mnt/ceph/jarredk/Cancer_Selection/GSE101961_normal_breast/GSE101961.unnorm.methyl.final.txt", quote=F, sep="\t")
write.table(mat.sig.pval.only, file="/mnt/ceph/jarredk/Cancer_Selection/GSE101961_normal_breast/GSE101961.signal.pvals.txt", quote=F, sep="\t")
#Normalize methylation data in same method as BSGS data

LT=function(X){
  
  Z=ifelse(is.na(X), NA, log(X/(1-X)))
  
  return(Z)
}

pseudocount=function(X){
  
  if(length(which(X==0))>0){ X[which(X==0)]=0.001 }
  if(length(which(X==1))>0){ X[which(X==1)]=0.999 }
  
  return(X)
  
}

tmethyl.M=apply(tmethyl.M, 2, pseudocount)
tmethyl.Mlt=apply(tmethyl.M, 2, LT)







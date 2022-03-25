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

#png("/mnt/ceph/jarredk/Cancer_Selection/TCGA.pca.plot.png")

#plot(TCGA.pca[,1], TCGA.pca[,2], color=as.factor(TCGA.meta$))


#####################################################################

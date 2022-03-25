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

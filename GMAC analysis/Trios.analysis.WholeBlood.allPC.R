library(stringi,lib="/mnt/lfs2/mdbadsha/Rpackages") #need for stri_sub ()
library(stringr,lib="/mnt/lfs2/mdbadsha/Rpackages")
library(MRPC,lib="/mnt/lfs2/mdbadsha/Rpackages")
library(WGCNA,lib="/mnt/lfs2/mdbadsha/Rpackages")
library(psych,lib="/mnt/lfs2/mdbadsha/Rpackages")
library(MatrixEQTL,lib="/mnt/lfs2/mdbadsha/Rpackages")

# cis-eQTL data from GTEx portal (https://gtexportal.org/)

WholeBlood.egenes.file.V8.1 <- read.delim (gzfile ("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt.gz"), header=TRUE)
dim(WholeBlood.egenes.file.V8.1)

# remove version number from gene id (e.g., ENSG00000227232.5 to ENSG00000227232)
WholeBlood.egenes.file.V8 <- WholeBlood.egenes.file.V8.1
WholeBlood.egenes.file.V8$gene_id <- gsub("\\..*","",WholeBlood.egenes.file.V8.1$gene_id)
WholeBlood.egenes.file.V8$gene_name <- gsub("\\..*","",WholeBlood.egenes.file.V8.1$gene_name)


# Read cis-eQTL file based on pvalue ascending data
# WholeBlood.egenes.file.accending.V8 <- WholeBlood.egenes.file.V8[order(as.numeric(as.character(WholeBlood.egenes.file.V8[,30]))) , ]
WholeBlood.egenes.file.accending.V8 <- WholeBlood.egenes.file.V8[order(WholeBlood.egenes.file.V8['pval_nominal']), ]
WholeBlood.all.Genes.V8 <- WholeBlood.egenes.file.accending.V8[,c(1,2,4,5,12,14,15,19,29,30)]
dim(WholeBlood.all.Genes.V8)
WholeBlood.all.Genes.V8[1:10,]


# To obtain the list of eGenes, select the rows with 'qval' ≤ 0.05.
# LC
WholeBlood.list.of.eGenes.V8 <- which(WholeBlood.all.Genes.V8[,9]<=0.05)
WholeBlood.eGenes.V8 <- WholeBlood.all.Genes.V8[WholeBlood.list.of.eGenes.V8,]
dim(WholeBlood.eGenes.V8)
WholeBlood.eGenes.V8[1:10,]

#unique cis-snp
WholeBlood.eGenes.V8.unique <- WholeBlood.eGenes.V8[!duplicated(WholeBlood.eGenes.V8$variant_id), ]
# cis-trans (LC)
#load(file="BrainCerebellar.eGenes.V8.RData")
LC.WholeBlood <- WholeBlood.eGenes.V8.unique[,c(5,1)]
colnames(LC.WholeBlood) <- c("snps", "cis.genes")
dim(LC.WholeBlood)
LC.WholeBlood[1:10,]


# load peer normalize data for WholeBlood 
load("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/PEER_Files_V8/Peerdata.WholeBlood.V8.RData")
dim(Peerdata.WholeBlood.V8)
Peerdata.WholeBlood.V8[1:10,1:5]

# rownames ENSG00000223972.5 to ENSG00000223972
Peerdata.WholeBlood.V8.1 <- Peerdata.WholeBlood.V8
# Peerdata.WholeBlood.V8.genes.names <- gsub("\\..*","",rownames(Peerdata.WholeBlood.V8)) 
rownames(Peerdata.WholeBlood.V8.1) <- gsub("\\..*","",rownames(Peerdata.WholeBlood.V8))
Peerdata.WholeBlood.V8.1[1:10,1:5]


# Gene expression data for trans-genes

# Read gene data from biomart
Genes.biomart.V8 <- read.delim2("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/Test_GTEx/mart_export.txt")
dim(Genes.biomart.V8)
Genes.biomart.V8[1:10,]

# Matching genes from biomart to GTEx v8
Match.gene.WholeBlood.V8 <- match(Genes.biomart.V8$Gene.stable.ID,rownames(Peerdata.WholeBlood.V8.1))
Match.gene.WholeBlood.V8.nona <- Match.gene.WholeBlood.V8[!is.na(Match.gene.WholeBlood.V8)] #deal with na
length(Match.gene.WholeBlood.V8.nona)

# Final gene expression data for trans-genes

WholeBlood.geneexpression.transdata.V8.1 <- Peerdata.WholeBlood.V8.1[Match.gene.WholeBlood.V8.nona,]
dim(WholeBlood.geneexpression.transdata.V8.1)
WholeBlood.geneexpression.transdata.V8.1[1:10,1:4]

save(WholeBlood.geneexpression.transdata.V8.1,file="/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/WholeBlood.geneexpression.transdata.V8.RData")

WholeBlood.geneexpression.transdata.V8 <- cbind(rownames(WholeBlood.geneexpression.transdata.V8.1), data.frame(WholeBlood.geneexpression.transdata.V8.1, row.names=NULL))
colnames(WholeBlood.geneexpression.transdata.V8)[1] <- "gene_id"
dim(WholeBlood.geneexpression.transdata.V8)
WholeBlood.geneexpression.transdata.V8[1:10,1:4]

# Save data for MatrixEQTL package
write.table(WholeBlood.geneexpression.transdata.V8, file = "/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/WholeBlood.geneexpression.transdata.V8.txt", sep = "\t",row.names = FALSE,col.names = TRUE)


# gene position data
Match.genepos.WholeBlood.V8 <- match(WholeBlood.geneexpression.transdata.V8$gene_id,Genes.biomart.V8$Gene.stable.ID)
Match.genepos.WholeBlood.V8.nona <- Match.genepos.WholeBlood.V8[!is.na(Match.genepos.WholeBlood.V8)] #deal with na
length(Match.genepos.WholeBlood.V8.nona)

# Final gene pos data for trans-genes
WholeBlood.genepos.V8 <- Genes.biomart.V8[Match.genepos.WholeBlood.V8.nona,c(1,5,3,4)]
dim(WholeBlood.genepos.V8)
colnames(WholeBlood.genepos.V8) <- c("gene_id" ,"chr","left","right")
WholeBlood.genepos.V8[1:10,]

# Save data for MatrixEQTL package
# save(WholeBlood.genepos.V8,file="/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/WholeBlood.genepos.V8.RData")

write.table(WholeBlood.genepos.V8, file = "/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/WholeBlood.genepos.V8.txt", sep = "\t",row.names = FALSE,col.names = TRUE)



# snp data processing for trans-genes
# load genotype data
load("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/snp.data.file.V8.RData")
dim(snp.data.file.V8)
snp.data.file.V8[1:10,1:8]
# match the snp from egene file with genotype file
match.snp.WholeBlood <- match(WholeBlood.eGenes.V8.unique$variant_id,snp.data.file.V8$ID)
match.snp.WholeBlood.nona <- match.snp.WholeBlood[!is.na(match.snp.WholeBlood)]
length(match.snp.WholeBlood.nona)

# match the sample from peer data file with genotype file
match.indv.WholeBlood <- match(colnames(Peerdata.WholeBlood.V8.1),colnames(snp.data.file.V8))
match.indv.WholeBlood.nona <- match.indv.WholeBlood[!is.na(match.indv.WholeBlood)]
length(match.indv.WholeBlood.nona)

WholeBlood.snp.data.V8.1 <- snp.data.file.V8[match.snp.WholeBlood.nona,match.indv.WholeBlood.nona]
rownames(WholeBlood.snp.data.V8.1) <- WholeBlood.eGenes.V8.unique$variant_id
WholeBlood.snp.data.V8.1[1:10,1:5]


# genotype data with first column is the snp id
WholeBlood.snp.data.V8 <- cbind(rownames(WholeBlood.snp.data.V8.1), data.frame(WholeBlood.snp.data.V8.1, row.names=NULL))
colnames(WholeBlood.snp.data.V8)[1] <- "snpid"
dim(WholeBlood.snp.data.V8)
WholeBlood.snp.data.V8[1:10,1:6]


# Save data for MatrixEQTL package
#save(WholeBlood.snp.data.V8.1, file = "/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/WholeBlood.snp.data.V8.RData")
write.table(WholeBlood.snp.data.V8, file = "/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/WholeBlood.snp.data.V8.txt", sep = "\t",row.names = FALSE,col.names = TRUE)



# snps pos data file (need to run MatrixEQTL package) snpid, chr and pos
WholeBlood.snpspos.V8.1 <- snp.data.file.V8[match.snp.WholeBlood,c(3,1,2)]
dim(WholeBlood.snpspos.V8.1)
colnames(WholeBlood.snpspos.V8.1) <- c("snpid", "chr","pos")
WholeBlood.snpspos.V8.1[1:10,]
# 

# convert chr1 to 1 (chr need same as genepos file)
WholeBlood.snpspos.V8 <- WholeBlood.snpspos.V8.1
WholeBlood.snpspos.V8$chr <- gsub('\\chr', '', WholeBlood.snpspos.V8.1$chr)
WholeBlood.snpspos.V8[1:10,]

# Save data for MatrixEQTL package
write.table(WholeBlood.snpspos.V8, file = "/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/WholeBlood.snpspos.V8.txt", sep = "\t",row.names = FALSE,col.names = TRUE)


# Run MatrixEQTL

useModel = modelLINEAR;
# Genotype file name
SNP_file_name = paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/WholeBlood.snp.data.V8.txt", sep="");
snps_location_file_name = paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/WholeBlood.snpspos.V8.txt", sep="");

# Gene expression file name
expression_file_name = paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/WholeBlood.geneexpression.transdata.V8.txt", sep="");
gene_location_file_name = paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/WholeBlood.genepos.V8.txt", sep="");

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 2e-2;
pvOutputThreshold_tra = 1e-5;


# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);


## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

Output.trans.cis.GTEx.WholeBlood.V8 <- Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);


save(Output.trans.cis.GTEx.WholeBlood.V8 ,file="/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/Output.trans.cis.GTEx.WholeBlood.V8.RData")

#Output.trans.cis.GTEx.WholeBlood.V8.useRData <- Output.trans.cis.GTEx.WholeBlood.V8

# Analysis step
# Total number of snp-trans
Output.trans.cis.GTEx.WholeBlood.V8$trans$neqtls

dim(Output.trans.cis.GTEx.WholeBlood.V8$trans$eqtls)


# unique trans-snp
length(unique(Output.trans.cis.GTEx.WholeBlood.V8$trans$eqtls[,1]))

# unique trans gene
length(unique(Output.trans.cis.GTEx.WholeBlood.V8$trans$eqtls[,2]))


dim(Output.trans.cis.GTEx.WholeBlood.V8$trans$eqtls)
LT.WholeBlood <- Output.trans.cis.GTEx.WholeBlood.V8$trans$eqtls[,c(1,2)]
colnames(LT.WholeBlood) <- c("snps", "trans.genes")
dim(LT.WholeBlood)
LT.WholeBlood[1:10,]


# Common snp between LC and LT
Common.snp.WholeBlood <- match(LC.WholeBlood$snps,LT.WholeBlood$snps)
Common.snp.WholeBlood.nona <- Common.snp.WholeBlood[!is.na(Common.snp.WholeBlood)]
length(Common.snp.WholeBlood.nona)


# Select TRANS GENES based on CIS genes
Common.snps.transgenes.WholeBlood <- LT.WholeBlood[Common.snp.WholeBlood.nona,]
dim(Common.snps.transgenes.WholeBlood)
Common.snps.transgenes.WholeBlood[1:5,]


# CIS GENES
Common.snp.cis.WholeBlood <- match(Common.snps.transgenes.WholeBlood$snps,LC.WholeBlood$snps)
Common.snp.cis.WholeBlood.nona <- Common.snp.cis.WholeBlood[!is.na(Common.snp.cis.WholeBlood)]
length(Common.snp.cis.WholeBlood.nona)

Common.snps.cisgenes.WholeBlood <- LC.WholeBlood[Common.snp.cis.WholeBlood.nona,]
dim(Common.snps.cisgenes.WholeBlood)
Common.snps.cisgenes.WholeBlood[1:5,]



# remove same cis and trans genes
Common.snps.cisgenes.WholeBlood.1 <- sapply(Common.snps.cisgenes.WholeBlood, as.character)
Common.snps.transgenes.WholeBlood.1 <- sapply(Common.snps.transgenes.WholeBlood, as.character)
same.LC.LT <- which(apply(Common.snps.cisgenes.WholeBlood.1==Common.snps.transgenes.WholeBlood.1,1,prod)==1)

if(length(same.LC.LT)!=0)
{
  Common.snps.cisgenes.WholeBlood <- Common.snps.cisgenes.WholeBlood.1[-c(same.LC.LT),]
  dim(Common.snps.cisgenes.WholeBlood)
  Common.snps.cisgenes.WholeBlood[1:5,]
  
  Common.snps.transgenes.WholeBlood <- Common.snps.transgenes.WholeBlood.1[-c(same.LC.LT),]
  dim(Common.snps.transgenes.WholeBlood)
  Common.snps.transgenes.WholeBlood[1:5,] 
}

if(length(same.LC.LT)==0){
  Common.snps.cisgenes.WholeBlood <- Common.snps.cisgenes.WholeBlood.1
  dim(Common.snps.cisgenes.WholeBlood)
  Common.snps.cisgenes.WholeBlood[1:5,]
  
  Common.snps.transgenes.WholeBlood <- Common.snps.transgenes.WholeBlood.1
  dim(Common.snps.transgenes.WholeBlood)
  Common.snps.transgenes.WholeBlood[1:5,] 
}


# formating snps-genes data

id <- dim(Common.snps.cisgenes.WholeBlood)[1] 

data.snp.cis.trans.WholeBlood <- NULL
cis.gene.WholeBlood <- NULL
trans.gene.WholeBlood <- NULL
snp.WholeBlood <- NULL

for (i in 1:id) {
  snp.WholeBlood[i] <- paste0(Common.snps.cisgenes.WholeBlood[i,1])
  cis.gene.WholeBlood[i] <- paste0(Common.snps.cisgenes.WholeBlood[i,2])
  trans.gene.WholeBlood[i] <- paste0(Common.snps.transgenes.WholeBlood[i,2])
  
  data.snp.cis.trans.WholeBlood[[i]]<-c(WholeBlood.snp.data.V8.1[which(rownames(WholeBlood.snp.data.V8.1)==snp.WholeBlood[i])[1],],
                                                            Peerdata.WholeBlood.V8.1[which(rownames(Peerdata.WholeBlood.V8.1)==cis.gene.WholeBlood[i])[1],],
                                                            Peerdata.WholeBlood.V8.1[which(rownames(Peerdata.WholeBlood.V8.1)==trans.gene.WholeBlood[i])[1],])
}

Col.names <- c(snp.WholeBlood, cis.gene.WholeBlood, trans.gene.WholeBlood)
data.snp.cis.trans.final.WholeBlood.V8.unique.snps <- matrix(as.numeric(unlist(data.snp.cis.trans.WholeBlood)),nrow=dim(Peerdata.WholeBlood.V8.1)[2],ncol=id*3, byrow = FALSE)

seq.list <- NULL
for(j in 1:id){
  seq.list[[j]]<-c(seq(j, id*3, id))
}
colnames(data.snp.cis.trans.final.WholeBlood.V8.unique.snps) <- Col.names[as.numeric(unlist(seq.list))]
rownames(data.snp.cis.trans.final.WholeBlood.V8.unique.snps) <- colnames(Peerdata.WholeBlood.V8.1)
dim(data.snp.cis.trans.final.WholeBlood.V8.unique.snps)
data.snp.cis.trans.final.WholeBlood.V8.unique.snps[1:10,1:6]

save(data.snp.cis.trans.final.WholeBlood.V8.unique.snps, file="/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/data.snp.cis.trans.final.WholeBlood.V8.unique.snps.RData")


Peerdata.WholeBlood2 <- t(Peerdata.WholeBlood.V8)
Peerdata.WholeBlood3 <- Peerdata.WholeBlood2[,apply(Peerdata.WholeBlood2, 2, var, na.rm=TRUE) != 0]
# pca matrix
PCs.WholeBlood <- prcomp(Peerdata.WholeBlood3,scale=TRUE)
PCs.matrix.WholeBlood <- PCs.WholeBlood$x

save(PCs.matrix.WholeBlood, file="/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/PCs.matrix.WholeBlood.RData")

data.all.unique.snps <- data.snp.cis.trans.final.WholeBlood.V8.unique.snps
#<- data.snp.cis.trans.final.WholeBlood.V8[,apply(data.snp.cis.trans.final.WholeBlood.V8, 2, var, na.rm=TRUE) != 0]
colnames(data.all.unique.snps) <- make.unique(colnames(data.all.unique.snps)) #duplicated colnames

All.Pvalues <- list()
List.significant.asso1.WholeBlood <- list()

# calculate associated PCs on all trios
for (m in 1:(dim(PCs.matrix.WholeBlood)[1]-1)) {
  corr.PCs <- corr.test(PCs.matrix.WholeBlood[,m],data.all.unique.snps,use = 'pairwise.complete.obs')
  # The p values
  Pvalues <- corr.PCs$p
  Pvalues.nona <- Pvalues[!is.na(Pvalues)]
  All.Pvalues [[m]] <- Pvalues.nona
  qobj <- qvalue(Pvalues.nona, fdr.level=0.10) 
  
  # Significant associations
  Significant.asso <- qobj$significant
  List.significant.asso1.WholeBlood[[m]] <- which(Significant.asso,useNames = TRUE)
}

save(List.significant.asso1.WholeBlood, file="/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/List.significant.asso1.WholeBlood.RData")

# apply MRPC with adjusted PCs on each trios

#truth for M0
#V1-->T1
Truth.M0 <- MRPCtruth$M0
Adj.M0<- as(Truth.M0,"matrix")
#V1-->T2
Adj.M01 <- matrix(0,nrow=3,ncol = 3)
rownames(Adj.M01) <- colnames(Adj.M01) <- colnames(Adj.M0)
Adj.M01[1,3] <- 1

#truth for M1
#V1-->T1-->T2
Truth.M1 <- MRPCtruth$M1
Adj.M1<- as(Truth.M1,"matrix")
#V1-->T2-->T1
Adj.M11 <- matrix(0,nrow=3,ncol = 3)
rownames(Adj.M11) <- colnames(Adj.M11) <- colnames(Adj.M1)
Adj.M11[1,3] <- 1
Adj.M11[3,2] <- 1

#truth for M2
#V1-->T1<--T2
Truth.M2 <- MRPCtruth$M2
Adj.M2<- as(Truth.M2,"matrix")
#V1-->T2<--T1
Adj.M21 <- matrix(0,nrow=3,ncol = 3)
rownames(Adj.M21) <-colnames(Adj.M21) <-colnames(Adj.M2)
Adj.M21[1,3] <- 1
Adj.M21[2,3] <- 1
#truth for M3
#V1-->T1, V1-->T2
Truth.M3 <- MRPCtruth$M3
Adj.M3 <- as(Truth.M3,"matrix")

#truth for M4
#V1-->T1, V1-->T2, T1--T2
Truth.M4 <- MRPCtruth$M4
Adj.M4 <- as(Truth.M4,"matrix")


#Initial results
M00.FDR.WholeBlood <- c()
M11.FDR.WholeBlood <- c()
M22.FDR.WholeBlood <- c()
M33.FDR.WholeBlood <- c()
M44.FDR.WholeBlood <- c()
List.Associated.PCs.WholeBlood <- list()

# columns number for trios
start.col <- seq(1,dim(data.snp.cis.trans.final.WholeBlood.V8.unique.snps)[2],3)
end.col <- seq(3,dim(data.snp.cis.trans.final.WholeBlood.V8.unique.snps)[2],3)

trios <- dim(data.snp.cis.trans.final.WholeBlood.V8.unique.snps)[2]/3
print(trios)

Match.PCs.col <- list()
List.Match.significant.trios <- list()

for (ii in 1:trios) {
  # data
  data <- data.snp.cis.trans.final.WholeBlood.V8.unique.snps[,start.col[ii]:end.col[ii]]
  # search in each PCs
  for (m1 in 1:(dim(PCs.matrix.WholeBlood)[1]-1)) {
    Match.PCs.col[[m1]] <- intersect(colnames(data),colnames(data.snp.cis.trans.final.WholeBlood.V8.unique.snps)[List.significant.asso1.WholeBlood[[m1]]])
  }
  List.Match.significant <- which(sapply(Match.PCs.col, function(e) length(e)!=0))
  List.Match.significant.trios[[ii]] <- List.Match.significant
  if(length(List.Match.significant)!=0)
  {
    data.withPC <- cbind(data,PCs.matrix.WholeBlood[,List.Match.significant])
    colnames(data.withPC)[4:ncol(data.withPC)] <- paste("PC",List.Match.significant,sep="")
  }
  else
  {
    data.withPC <- data 
  }
  
  
  
  n <- nrow (data.withPC)
  V <- colnames(data.withPC)     # Column names
  
  # Classical correlation
  suffStat <- list(C = cor(data.withPC,use = "complete.obs"),
                   n = n)
  
  MRPC.fit.FDR <- MRPC(data.withPC,
                       suffStat,
                       GV = 1,
                       FDR = 0.05,
                       indepTest = 'gaussCItest',
                       labels = V,
                       FDRcontrol = TRUE,
                       verbose = FALSE)
  
  #plot(MRPC.fit_FDR)
  Adj.infe1 <- as( MRPC.fit.FDR@graph,"matrix")
  Adj.infe <- Adj.infe1[1:3,1:3] #only consider snp, cis, trans
  colnames(Adj.infe) <- rownames(Adj.infe) <- colnames(Adj.M01) 
  
  if(identical(Adj.M0,Adj.infe) || identical(Adj.M01,Adj.infe)){
    M00.FDR.WholeBlood[ii]<-ii
  }
  
  if(identical(Adj.M1,Adj.infe) || identical(Adj.M11,Adj.infe)){
    M11.FDR.WholeBlood[ii]<-ii
  }
  if(identical(Adj.M2,Adj.infe) || identical(Adj.M21,Adj.infe)){
    M22.FDR.WholeBlood[ii]<-ii
  }
  if(identical(Adj.M3,Adj.infe)){
    M33.FDR.WholeBlood[ii]<-ii
  }
  if(identical(Adj.M4,Adj.infe)){
    M44.FDR.WholeBlood[ii]<-ii
  }
  
}

save(List.Match.significant.trios,file="/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/List.Match.significant.trios.RData")



List.models.WholeBlood.all <- list()
List.models.WholeBlood.all$M0 <- M00.FDR.WholeBlood[!is.na(M00.FDR.WholeBlood)] 
List.models.WholeBlood.all$M1 <- M11.FDR.WholeBlood[!is.na(M11.FDR.WholeBlood)] 
List.models.WholeBlood.all$M2 <- M22.FDR.WholeBlood[!is.na(M22.FDR.WholeBlood)] 
List.models.WholeBlood.all$M3 <- M33.FDR.WholeBlood[!is.na(M33.FDR.WholeBlood)] 
List.models.WholeBlood.all$M4 <- M44.FDR.WholeBlood[!is.na(M44.FDR.WholeBlood)] 

save(List.models.WholeBlood.all,file="/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/List.models.WholeBlood.all.RData")

All.models.WholeBlood.all<- matrix(c(length(List.models.WholeBlood.all$M0),
                                                         length(List.models.WholeBlood.all$M1),
                                                         length(List.models.WholeBlood.all$M2),
                                                         length(List.models.WholeBlood.all$M3),
                                                         length(List.models.WholeBlood.all$M4),
                                                         length(setdiff(1:trios,c(List.models.WholeBlood.all$M0,List.models.WholeBlood.all$M1,
                                                                                  List.models.WholeBlood.all$M2,List.models.WholeBlood.all$M3,
                                                                                  List.models.WholeBlood.all$M4))),
                                                         length(WholeBlood.eGenes.V8.unique$variant_id),
                                                         dim(Peerdata.WholeBlood.V8)[1],
                                                         trios,dim(Peerdata.WholeBlood.V8)[2],
                                                         dim(PCs.matrix.WholeBlood)[2]-1),
                                                       nrow = 1,ncol = 11)
colnames(All.models.WholeBlood.all) <- c("M0","M1","M2","M3","M4","Others","snps","genes","trios","samples","no.PCs")

print(All.models.WholeBlood.all)

save(All.models.WholeBlood.all ,file="/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/All.models.WholeBlood.all.RData")


#List.Match.significant.trios.1 <- which(sapply(List.Match.significant.trios, function(e) length(e)!=0))
PCasso.pc.all.trios <- which(sapply(List.Match.significant.trios, function(e) length(e)!=0))
length(PCasso.pc.all.trios)
#summary of associated PCs
Match.significant.trios.summary <- lengths(List.Match.significant.trios, use.names = TRUE)
summary(Match.significant.trios.summary)
#Total no of pc associated
List.significant.asso2 <- which(sapply(List.significant.asso1.WholeBlood, function(e) length(e)!=0))
length(List.significant.asso2)

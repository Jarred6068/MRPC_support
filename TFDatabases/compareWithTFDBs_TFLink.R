##################################################
# 1. Load functions and set up the environment
##################################################
source("functions_TFDatabases.R")
library(openxlsx)
library(data.table)

#################################################################
# 2. Read in mediation trios
# Combine all the cis-gene (or trans-gene) mediation trios across
# tissues into one data frame
#################################################################
setwd("~/Documents/CisTrans/Data/SuppTables_202412")
path <- "Supplemental_Table_S8_MRPC_LOND_M1_cis_each_tissue.xlsx"
trios.all.m11 <- readMRGNTrios (path)
dim (trios.all.m11)
# 395  26

path <- "Supplemental_Table_S9_MRPC_LOND_M1_trans_each_tissue.xlsx"
trios.all.m12 <- readMRGNTrios (path)
dim (trios.all.m12)
# 1323   26

##################################################
# 3. Extract TF-target pairs from final trios with 
# matched TFs
# TFLink downloaded here:
# https://tflink.net/download/
# under 'Homo sapiens small and large-scale interaction table'
# filename: TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv.gz
##################################################
tflink <- read.table(gzfile("TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv.gz"), header = TRUE, sep = "\t")
dim (tflink)
# 6739357      15
all_unique <- unique(tflink$Name.TF)
length (all_unique)
# 1606


# Extract overlapping TFs from TFLink
# Extract the corresponding trios
# m11
tmp <- match (all_unique, trios.all.m11$Cis.Gene.Name)
matched.tfs.ids <- which (!is.na(tmp))
tf.name.m11 <- all_unique[matched.tfs.ids]
# find trios with TF
# using TFLink
tf.targets.tflink.m11 <- NULL
for (i in 1:length(tf.name.m11)) {
  print (i)
  tf.targets.tflink.m11[[i]] <- tflink$Name.Target[which ((tflink$Name.TF==tf.name.m11[i]))]
}
result.m11.tflink <- extractTFsFromTrios(trios.all.m11, tf.name.m11, is.cis=TRUE)
write.table (result.m11.tflink, "trios_m11_cisgene_as_tf_tflink.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

unlist (lapply (tf.targets.tflink.m11, length))
#[1]  145    2 6012    4    1    1 2205  439 1372  410  682  835  764
sum (unlist (lapply (tf.targets.tflink.m11, length)))
#[1] 12872

# m12
tmp <- match (all_unique, trios.all.m12$Trans.Gene.Name)
matched.tfs.ids <- which (!is.na(tmp))
tf.name.m12 <- all_unique[matched.tfs.ids]
# find trios with TF
# using TFLink
tf.targets.tflink.m12 <- NULL
for (i in 1:length(tf.name.m12)) {
  print (i)
  tf.targets.tflink.m12[[i]] <- tflink$Name.Target[which ((tflink$Name.TF==tf.name.m12[i]))]
}
result.m12.tflink <- extractTFsFromTrios(trios.all.m12, tf.name.m12, is.cis=FALSE)
write.table (result.m12.tflink, "trios_m12_transgene_as_tf_tflink.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)


##################################################
# 4. identify candidate trios of the matched TF-target pairs
##################################################
# Set up Entrez and Ensembl ID databases
library(org.Hs.eg.db, lib='/mnt/ceph/audreyf/Rpackages')
db.entrez <- as.list(org.Hs.egALIAS2EG)

x <- org.Hs.egENSEMBL
mapped_genes <- mappedkeys(x)
db.ensembl <- as.list(x[mapped_genes])


# Read in the candidate trios
tissues.vec <- read.csv("tissuenames.csv", header = T)
result.m11.tflink.trios <- extractTFTargetFromGTEx(tf.name.m11, tf.targets.tflink.m11, tissues=tissues.vec$tissue.name1, file.path='/mnt/zeta/audreyf/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/', db.entrez, db.ensembl)
write.table (result.m11.tflink.trios, "result_m11_matched_candidate_trios_tflink.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)

result.m12.tflink.trios <- extractTFTargetFromGTEx(tf.name.m12, tf.targets.tflink.m12, tissues=tissues.vec$tissue.name1, file.path='/mnt/zeta/audreyf/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/', db.entrez, db.ensembl)
write.table (result.m12.tflink.trios, "result_m12_matched_candidate_trios_tflink.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)


##################################################
# 5. Find same TF-target pairs in GTEx trios and the TF database
##################################################
compareTables (result.m11.tflink, result.m11.tflink.trios)
compareTables (result.m12.tflink, result.m12.tflink.trios, is.cis = FALSE)

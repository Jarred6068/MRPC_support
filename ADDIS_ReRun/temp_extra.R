

n <- nrow (trio.test1)
V <- colnames(trio.test1)     # Column names

# Classical correlation
suffStat <- list(C = cor(trio.test1,use = "complete.obs"),
                 n = n)

MRPC.fit.FDR.addis <- MRPC(trio.test1,
                     suffStat,
                     GV = 1,
                     FDR = 0.05,
                     indepTest = 'gaussCItest',
                     labels = V,
                     FDRcontrol = "ADDIS",
                     verbose = FALSE)


MRPC.fit.FDR.lond <- MRPC(trio.test1,
                     suffStat,
                     GV = 1,
                     FDR = 0.05,
                     indepTest = 'gaussCItest',
                     labels = V,
                     FDRcontrol = "LOND",
                     verbose = FALSE)

print(as(MRPC.fit.FDR.addis@graph, "matrix"))
print(as(MRPC.fit.FDR.lond@graph, "matrix"))


#snp.data[, snp*3+(1:3)]

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


#read in the data
snps=loadRData(file="/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/data.snp.cis.trans.final.WholeBlood.V8.unique.snps.RData")
pc.matrix=loadRData("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/PCs.matrix.WholeBlood.RData")
sig.pc=loadRData("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/List.Match.significant.trios.RData")
list.LOND=loadRData("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/List.models.WholeBlood.all.RData")
list.ADDIS=loadRData("/mnt/ceph/jarredk/AddisReRunFiles/List.models.WholeBlood.all.RData")

#WholeBlood Example
list.LOND$M1
length(list.LOND$M1)

list.ADDIS$M1
length(list.ADDIS$M1)

#see which elements in lond$M1 are also in addis$M1

mat1=match(list.LOND$M1, list.ADDIS$M1)

#the elements returned as NA are found in LOND$M1 but not ADDIS$M1

#return the matched and unmatched elements
matched.elements=list.LOND$M1[!is.na(mat1)]
unmatched.elements=list.LOND$M1[is.na(mat1)]
print(matched.elements, unmatched.elements)










#==============================================================================================
#-----------------------------------Example--Code--use-----------------------------------------
#==============================================================================================
source("/mnt/ceph/jarredk/ADDIS_verify/ADDIS_Post_Analysis_processing.R")
#for reading in the tissue names
tissue.names=read.csv("/mnt/ceph/jarredk/AddisReRunFiles/tissuenames.csv", header = T)

#catalog changes from LOND-->ADDIS for tissue "Whole Blood" and "Prostate"
##note: tissue names follow conventions tissue.names[,2] i.e those for the _AllPC/ 
##directory path names. Therefore the second column of tissuenames.csv can be
##indexed and used directly as input

out1=ADDIS.PostProc(tissue.names = c("WholeBlood", "Prostate"))
##the output has 3 fields:
summary(out1)
#--$tiss--contains the catalog model/class changes for all desired tissue types  
summary(out1$tiss)
summary(out1$tiss$WholeBlood)
summary(out1$tiss$WholeBlood$Lond.M1)
#--$diags--the trios classified as the same model between FDR methods
summary(out1$diags)
#--$all.tables--A table of the the number of trios classified differently and the
##              location (as model type) under ADDIS. the diagonals of the table 
##              are the diagonal values (congruent trio class.) from $diags
summary(out1$all.tables)
#check the table for WholeBlood
print(out1$all.tables$WholeBlood)

#the function Lond2Addis.lookup() provides a convienent way to look deeper at specific
#trios whose model designation changed from LOND-->ADDIS

#picking a random trio in WholeBlood
out1$tiss$WholeBlood$Lond.M1

#lets pick the first element that was classified as M4
out2=Lond2Addis.lookup(trio.index=out1$tiss$WholeBlood$Lond.M1$Addis.M4[1], tissue.name = "WholeBlood")

#NOTE: the option r.pack.lib=  gives the path name for you Rpackages (or where ever you have MRPC). you 
#must specify this options with a character string containing the path name for the package MRPC

#In the most basic sense with all options set ==TRUE, the function outputs a list with 4 fields. The first 
#is the matrix of data for the trio (with or without pc's depending on options).
summary(out2)
print(head(out2$trio.with.pc))
#The second element is the correlation matrix for the trio. 
print(out2$correlation)
#And the 3rd and 4th are the matrix of the graph from MRPC run using LOND and ADDIS on the trio
print(out2$MRPC.table.LOND)
print(out2$MRPC.table.ADDIS)










#==========================================================================================
#-----------------------------------CHECKING-M1's------------------------------------------
#==========================================================================================

source("/mnt/ceph/jarredk/ADDIS_verify/ADDIS_Post_Analysis_processing.R")

#check data for Uterus, WholeBlood,Thyroid
out1.1=ADDIS.M1.check(tissue.names=c("Uterus","WholeBlood","Thyroid"))

#-------------------------------Uterus-------------------------------------
#M1 class that stayed M1 in Addis 
matched.1 = match(out1.1$Catalog$tiss$Uterus$ADDIS$type1, out1.1$Catalog$tiss$Uterus$LOND$type1)
matched.2 = match(out1.1$Catalog$tiss$Uterus$ADDIS$type2, out1.1$Catalog$tiss$Uterus$LOND$type2)
#find the trios that did not change 
incommon.M1.1=out1.1$Catalog$tiss$Uterus$ADDIS$type1[!is.na(matched.1)]
incommon.M1.2=out1.1$Catalog$tiss$Uterus$ADDIS$type2[!is.na(matched.2)]
#display
print(stayed.M1.1)
print(stayed.M1.2)
#find the trios that are unique to ADDIS for Type 1 $ 2
M1.A.unique.type1=out1.1$Catalog$tiss$Uterus$ADDIS$type1[is.na(matched.1)]
M1.A.unique.type2=out1.1$Catalog$tiss$Uterus$ADDIS$type2[is.na(matched.2)]
#display
print(M1.A.unique.type1)
print(M1.A.unique.type2)

#find the trios that are unique to LOND for Type 1 $ 2
M1.L.unique.type1=out1.1$Catalog$tiss$Uterus$LOND$type1[-c(na.omit(matched.1))]
M1.L.unique.type2=out1.1$Catalog$tiss$Uterus$LOND$type2[-c(na.omit(matched.2))]
#display
print(M1.L.unique.type1)
print(M1.L.unique.type2)


catdata=ADDIS.PostProc(tissue.names=c("Uterus","WholeBlood","Thyroid"))

catdata$tiss$Uterus$Lond.M1
M1.L.unique.type1
M1.L.unique.type2

list.LOND.Uterus=loadRData("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/Uterus_AllPC/List.models.Uterus.all.RData")
list.ADDIS.Uterus=loadRData("/mnt/ceph/jarredk/AddisReRunFiles/List.models.Uterus.all.RData")

# location.ADDIS=list()
# for(i in 1:length(list.LOND.Uterus)){
#   
#   indicies=match(M1.L.unique.type2, list.ADDIS.Uterus[[i]])
#   indicies=indicies[!is.na(indicies)]
#   print(indicies)
#   location.ADDIS[[i]]=list.ADDIS.Uteus[[i]][indicies]
#   
# }


lookup1=Lond2Addis.lookup(trio.index=M1.L.unique.type2[1], tissue.name=c("Uterus"))
lookup1[2:4]
trio=lookup1$trio.with.pc
summary(lm(trio[,3]~trio[,-3]))
summary(lm(trio[,2]~trio[,-2]))


lookup2=Lond2Addis.lookup(trio.index=incommon.M1.1[1], tissue.name=c("Uterus"))
lookup2[2:4]
trio=lookup2$trio.with.pc
summary(lm(trio[,3]~trio[,-3]))
summary(lm(trio[,2]~trio[,-2]))









#HIC and resampling:

source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")

#check all M1's for CellsEBVtransformedlymphocytes
M1.lymphocytes=ADDIS.M1.check("CellsEBVtransformedlymphocytes")

indexes2=c(M1.lymphocytes$Catalog$tiss$CellsEBVtransformedlymphocytes$ADDIS$type2)

#M1.lymphocytes$Catalog$tiss$CellsEBVtransformedlymphocytes$ADDIS$type1

triosM1.hic.result=interaction_check(hic.filename="/mnt/ceph/jarredk/HiC_Analyses/lymphoblastoid_cells/ENCFF355OWW.hic", 
                                     trios=indexes2,
                                     resolution=10000,
                                     search.size=100000,
                                     tiss="CellsEBVtransformedlymphocytes")
data1=na.omit(triosM1.hic.result)
str(data1)
output=Resample_interactions(filePath = "/mnt/ceph/jarredk/HiC_Analyses/fibroblast_cells/ENCFF768UBD.hic",
                      chrs=c(paste(data1$cis.chr[2]),paste(data1$trans.chr[2])),
                      res=10000,
                      search.size=100000,
                      plot.hist=TRUE)


output=Resample_interactions(filePath = "/mnt/ceph/jarredk/HiC_Analyses/fibroblast_cells/ENCFF768UBD.hic",
                             chrs=c(paste(data1$cis.chr[6]),paste(data1$trans.chr[6])),
                             res=10000,
                             search.size=100000)


save(triosM1.hic.result, file = "/mnt/ceph/jarredk/HiC_Analyses/lymphoblastoid_cells/resampled.lymph.Rdata")
save(triosM1.hic.result2, file = "/mnt/ceph/jarredk/HiC_Analyses/Lung/resampled.Lung.Rdata")
save(triosM1.hic.result3, file = "/mnt/ceph/jarredk/HiC_Analyses/Skin/resampled.Skin.Rdata")
save(triosM1.hic.result4, file = "/mnt/ceph/jarredk/HiC_Analyses/fibroblast_cells/resampled.fibro.Rdata")


install=FALSE
if(install==TRUE){
#two current HiC data sets are:
#(1). GSE128678_processed_files.tar.gz ---> Lymphoblastoid Cell Lines via GEO Repo
    # related file: /mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/CellsEBVtransformedlymphocytes_AllPC/CellsEBVtransformedlymphocytes.eGenes.V8.unique.RData

#(2). ENCFF355OWW.hic ---> Homo sapiens: GM12878

    #related files

#install HiCcompare package from BiocManager:
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager", lib = "/mnt/ceph/jarredk/Rpackages")
BiocManager::install("HiCcompare", lib="/mnt/ceph/jarredk/Rpackages")
BiocManager::install("HiTC", lib="/mnt/ceph/jarredk/Rpackages")
BiocManager::install("HiCDataHumanIMR90", lib="/mnt/ceph/jarredk/Rpackages")
library(HiCcompare)
browseVignettes("HiTC")
#Straw is an API for fast data extraction for .hic files that provides programmatic access to the matrices.
#Straw is compatible with R via the Rcpp library

library(Rcpp)
sourceCpp("straw-R.cpp")

#E.G:

#      A<-straw_R("NONE drosophila.hic arm_2L arm_2L BP 100000")

#In the above example, A is a data frame containing the counts information:
#To extract the raw matrix corresponding to chromosome 22 at 500kb resolution 
#we would use the following command within the terminal:

#      straw_R("NONE GSE63525_K562_combined_30.hic  22 22 BP 500000")


#we can pair specific chromosome location and resolutions with data contained in the snp data in
#/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC




#method 2 of using straw !!this method verified to work!!
install.packages('remotes', lib="/mnt/ceph/jarredk/Rpackages")

remotes::install_github("aidenlab/straw/R")
hic.data.frame <- strawr::straw("KR", "/path/to/file.hic", "11", "11", "BP", 10000)

}



















#ensembl.org


#=============================================================================================
#-----------------------------------------Example_HiC_Code------------------------------------
#=============================================================================================

source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")
tissue.names=read.csv("/mnt/ceph/jarredk/AddisReRunFiles/tissuenames.csv", header = T)
print(tissue.names[22,2])
M1.lymphocytes=ADDIS.M1.check(tissue.names[22,2])
indexes=M1.lymphocytes$Catalog$tiss$CellsEBVtransformedlymphocytes$ADDIS$type2[1:4]

trio.attributes=get_trio_attr(trio.index=indexes, tissue.name=tissue.names[22,2], verbose=TRUE)
str(trio.attributes)

trio.attributes$Attributes$trans

trio.attributes$Attributes$cis

hic1=extract_hic(fileName="/mnt/ceph/jarredk/HiC_Analyses/lymphoblastoid_cells/ENCFF355OWW.hic", 
                 chrs=c("22","11"), 
                 resol=5000)

min(hic1$x)
max(hic1$x)
min(hic1$y)
max(hic1$y)



hic2=extract_hic(fileName="/mnt/ceph/jarredk/HiC_Analyses/lymphoblastoid_cells/ENCFF355OWW.hic", 
                 chrs=c("2","19"), 
                 resol=5000)

min(hic2$x)
max(hic2$x)
min(hic2$y)
max(hic2$y)

#example using interaction_check
interaction_check(hic.filename="/mnt/ceph/jarredk/HiC_Analyses/lymphoblastoid_cells/ENCFF355OWW.hic", 
                  trios=indexes,
                  resolution=10000,
                  search.size=100000,
                  tiss="CellsEBVtransformedlymphocytes")

#=====================================================================================================
#--------------------------------------END_Example_Codes----------------------------------------------
#=====================================================================================================





#=====================================================================================================
#------------------------------Run-Analysis-For-Available-Tissues-------------------------------------
#=====================================================================================================
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


#check all M1's for Lung
M1.Lung=ADDIS.M1.check("Lung")

indexes3=c(M1.Lung$Catalog$tiss$Lung$ADDIS$type2)

triosM1.hic.result2=interaction_check(hic.filename="/mnt/ceph/jarredk/HiC_Analyses/Lung/ENCFF366ERB.hic", 
                                     trios=indexes3,
                                     resolution=10000,
                                     search.size=100000,
                                     tiss="Lung")


#check all M1's skin
M1.skin=ADDIS.M1.check("SkinNotSunExposed")

indexes3=c(M1.fibroblasts$Catalog$tiss$SkinNotSunExposed$ADDIS$type2)

triosM1.hic.result3=interaction_check(hic.filename="/mnt/ceph/jarredk/HiC_Analyses/Skin/ENCFF569RJM.hic", 
                                     trios=indexes3,
                                     resolution=10000,
                                     search.size=50000,
                                     tiss="SkinNotSunExposed")


#check all M1's for Fibroblasts
M1.fibroblasts=ADDIS.M1.check("CellsCulturedfibroblasts")

indexes3=c(M1.fibroblasts$Catalog$tiss$CellsCulturedfibroblasts$ADDIS$type2)

triosM1.hic.result4=interaction_check(hic.filename="/mnt/ceph/jarredk/HiC_Analyses/fibroblast_cells/ENCFF768UBD.hic", 
                                     trios=indexes3,
                                     resolution=10000,
                                     search.size=100000,
                                     tiss="CellsCulturedfibroblasts")

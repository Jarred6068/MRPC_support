
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
#---------------------------Run-Analysis-For-Available-Tissues::ADDIS---------------------------------
#=====================================================================================================
#set.seed(566)

source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")

#check all M1's for CellsEBVtransformedlymphocytes
M1.lymphocytes=ADDIS.M1.check("CellsEBVtransformedlymphocytes")

indexes=c(M1.lymphocytes$Catalog$tiss$CellsEBVtransformedlymphocytes$ADDIS$type2)

#M1.lymphocytes$Catalog$tiss$CellsEBVtransformedlymphocytes$ADDIS$type1

triosM1.hic.result=interaction_check(hic.filename="/mnt/ceph/jarredk/HiC_Analyses/Data/lymphoblastoid_cells/ENCFF355OWW.hic", 
                                     trios=indexes,
                                     resolution=10000,
                                     search.size=100000,
                                     tiss="CellsEBVtransformedlymphocytes",
                                     verbose=TRUE,
                                     plot.h=TRUE,
                                     FDR="ADDIS")

write.csv(triosM1.hic.result$summary.table1, file = "/mnt/ceph/jarredk/HiC_Analyses/ADDIS/lymphoblastoid_cells/HiC_lymphocyte_result_ADDIS.csv")
save(triosM1.hic.result, file = "/mnt/ceph/jarredk/HiC_Analyses/ADDIS/lymphoblastoid_cells/resampled.lymph.Rdata")
triosM1.hic.result$summary.table1


#check all M1's for Lung
M1.Lung=ADDIS.M1.check("Lung")

indexes2=c(M1.Lung$Catalog$tiss$Lung$ADDIS$type2)

triosM1.hic.result2=interaction_check(hic.filename="/mnt/ceph/jarredk/HiC_Analyses/Data/Lung/ENCFF366ERB.hic", 
                                     trios=indexes2,
                                     resolution=10000,
                                     search.size=100000,
                                     tiss="Lung",
                                     verbose=TRUE,
                                     plot.h=TRUE,
                                     FDR="ADDIS")

write.csv(triosM1.hic.result2$summary.table1, file = "/mnt/ceph/jarredk/HiC_Analyses/ADDIS/Lung/HiC_Lung_result_ADDIS.csv")
save(triosM1.hic.result2, file = "/mnt/ceph/jarredk/HiC_Analyses/ADDIS/Lung/resampled.Lung.Rdata")
triosM1.hic.result2$summary.table1

#check all M1's skin
M1.skin=ADDIS.M1.check("SkinNotSunExposed")

indexes3=c(M1.skin$Catalog$tiss$SkinNotSunExposed$ADDIS$type2)

triosM1.hic.result3=interaction_check(hic.filename="/mnt/ceph/jarredk/HiC_Analyses/Data/Skin/ENCFF569RJM.hic", 
                                     trios=indexes3,
                                     resolution=10000,
                                     search.size=100000,
                                     tiss="SkinNotSunExposed",
                                     verbose=TRUE,
                                     plot.h=TRUE,
                                     FDR="ADDIS")

write.csv(triosM1.hic.result3$summary.table1, file = "/mnt/ceph/jarredk/HiC_Analyses/ADDIS/Skin/HiC_Skin_result_ADDIS.csv")
save(triosM1.hic.result3, file = "/mnt/ceph/jarredk/HiC_Analyses/ADDIS/Skin/resampled.Skin.Rdata")
triosM1.hic.result3$summary.table1

#check all M1's for Fibroblasts
M1.fibroblasts=ADDIS.M1.check("CellsCulturedfibroblasts")

indexes4=c(M1.fibroblasts$Catalog$tiss$CellsCulturedfibroblasts$ADDIS$type2)

triosM1.hic.result4=interaction_check(hic.filename="/mnt/ceph/jarredk/HiC_Analyses/Data/fibroblast_cells/ENCFF768UBD.hic", 
                                     trios=indexes4,
                                     resolution=10000,
                                     search.size=100000,
                                     tiss="CellsCulturedfibroblasts",
                                     verbose=TRUE,
                                     plot.h=TRUE,
                                     FDR="ADDIS")


write.csv(triosM1.hic.result4$summary.table1, file = "/mnt/ceph/jarredk/HiC_Analyses/ADDIS/fibroblast_cells/HiC_Fibroblasts_result_ADDIS.csv")
save(triosM1.hic.result4, file = "/mnt/ceph/jarredk/HiC_Analyses/ADDIS/fibroblast_cells/resampled.fibro.ADDIS.Rdata")
triosM1.hic.result4$summary.table1






#=====================================================================================================
#---------------------------Run-Analysis-For-Available-Tissues::LOND----------------------------------
#=====================================================================================================
source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")

#check all M1's for CellsEBVtransformedlymphocytes
M1.lymphocytes=ADDIS.M1.check("CellsEBVtransformedlymphocytes")

indexes=c(M1.lymphocytes$Catalog$tiss$CellsEBVtransformedlymphocytes$LOND$type2)

#M1.lymphocytes$Catalog$tiss$CellsEBVtransformedlymphocytes$ADDIS$type1

triosM1.hic.result=interaction_check(hic.filename="/mnt/ceph/jarredk/HiC_Analyses/Data/lymphoblastoid_cells/ENCFF355OWW.hic", 
                                     trios=indexes,
                                     resolution=10000,
                                     search.size=100000,
                                     tiss="CellsEBVtransformedlymphocytes",
                                     verbose=TRUE,
                                     plot.h=TRUE,
                                     FDR="LOND")

write.csv(triosM1.hic.result$summary.table1, file = "/mnt/ceph/jarredk/HiC_Analyses/LOND/lymphoblastoid_cells/HiC_lymphocyte_result_LOND.csv")
save(triosM1.hic.result, file = "/mnt/ceph/jarredk/HiC_Analyses/LOND/lymphoblastoid_cells/resampled.lymph.LOND.Rdata")
triosM1.hic.result$summary.table1
triosM1.hic.result$summary.table2


#check all M1's for Lung
M1.Lung=ADDIS.M1.check("Lung")

indexes2=c(M1.Lung$Catalog$tiss$Lung$LOND$type2)

triosM1.hic.result2=interaction_check(hic.filename="/mnt/ceph/jarredk/HiC_Analyses/Data/Lung/ENCFF366ERB.hic", 
                                      trios=indexes2,
                                      resolution=10000,
                                      search.size=100000,
                                      tiss="Lung",
                                      verbose=TRUE,
                                      plot.h=TRUE,
                                      FDR="LOND")

write.csv(triosM1.hic.result2$summary.table1, file = "/mnt/ceph/jarredk/HiC_Analyses/LOND/Lung/HiC_Lung_result_LOND.csv")
save(triosM1.hic.result2, file = "/mnt/ceph/jarredk/HiC_Analyses/LOND/Lung/resampled.Lung.LOND.Rdata")
triosM1.hic.result2$summary.table1

#check all M1's skin
M1.skin=ADDIS.M1.check("SkinNotSunExposed")

indexes3=c(M1.skin$Catalog$tiss$SkinNotSunExposed$LOND$type2)

triosM1.hic.result3=interaction_check(hic.filename="/mnt/ceph/jarredk/HiC_Analyses/Data/Skin/ENCFF569RJM.hic", 
                                      trios=indexes3,
                                      resolution=10000,
                                      search.size=100000,
                                      tiss="SkinNotSunExposed",
                                      verbose=TRUE,
                                      plot.h=TRUE,
                                      FDR="LOND")

write.csv(triosM1.hic.result3$summary.table1, file = "/mnt/ceph/jarredk/HiC_Analyses/LOND/Skin/HiC_Skin_result_LOND.csv")
save(triosM1.hic.result3, file = "/mnt/ceph/jarredk/HiC_Analyses/LOND/Skin/resampled.Skin.LOND.Rdata")
triosM1.hic.result3$summary.table1

#check all M1's for Fibroblasts
M1.fibroblasts=ADDIS.M1.check("CellsCulturedfibroblasts")

indexes4=c(M1.fibroblasts$Catalog$tiss$CellsCulturedfibroblasts$LOND$type2)

triosM1.hic.result4=interaction_check(hic.filename="/mnt/ceph/jarredk/HiC_Analyses/Data/fibroblast_cells/ENCFF768UBD.hic", 
                                      trios=indexes4,
                                      resolution=10000,
                                      search.size=100000,
                                      tiss="CellsCulturedfibroblasts",
                                      verbose=TRUE,
                                      plot.h=TRUE,
                                      FDR="LOND")


write.csv(triosM1.hic.result4$summary.table1, file = "/mnt/ceph/jarredk/HiC_Analyses/LOND/fibroblast_cells/HiC_Fibroblasts_result_LOND.csv")
save(triosM1.hic.result4, file = "/mnt/ceph/jarredk/HiC_Analyses/LOND/fibroblast_cells/resampled.fibro.LOND.Rdata")
triosM1.hic.result4$summary.table1



                      







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

write.csv(triosM1.hic.result$summary.table1, file = "/mnt/ceph/jarredk/HiC_Analyses/ADDIS/lymphoblastoid_cells/HiC_lymphocyte_result_ADDIS.csv",
          row.names = F)
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

write.csv(triosM1.hic.result2$summary.table1, file = "/mnt/ceph/jarredk/HiC_Analyses/ADDIS/Lung/HiC_Lung_result_ADDIS.csv",
          row.names = F)
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

write.csv(triosM1.hic.result3$summary.table1, file = "/mnt/ceph/jarredk/HiC_Analyses/ADDIS/Skin/HiC_Skin_result_ADDIS.csv",
          row.names = F)
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


write.csv(triosM1.hic.result4$summary.table1, file = "/mnt/ceph/jarredk/HiC_Analyses/ADDIS/fibroblast_cells/HiC_Fibroblasts_result_ADDIS.csv",
          row.names = F)
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

write.csv(triosM1.hic.result$summary.table1, file = "/mnt/ceph/jarredk/HiC_Analyses/LOND/lymphoblastoid_cells/HiC_lymphocyte_result_LOND.csv",
          row.names = F)
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

write.csv(triosM1.hic.result2$summary.table1, file = "/mnt/ceph/jarredk/HiC_Analyses/LOND/Lung/HiC_Lung_result_LOND.csv", 
          row.names = F)
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

write.csv(triosM1.hic.result3$summary.table1, file = "/mnt/ceph/jarredk/HiC_Analyses/LOND/Skin/HiC_Skin_result_LOND.csv",
          row.names = F)
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


write.csv(triosM1.hic.result4$summary.table1, file = "/mnt/ceph/jarredk/HiC_Analyses/LOND/fibroblast_cells/HiC_Fibroblasts_result_LOND.csv",
          row.names = F)
save(triosM1.hic.result4, file = "/mnt/ceph/jarredk/HiC_Analyses/LOND/fibroblast_cells/resampled.fibro.LOND.Rdata")
triosM1.hic.result4$summary.table1





source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")


get_trio_attr(indexes[1:10])


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


triosM1.hic.result$summary.table1

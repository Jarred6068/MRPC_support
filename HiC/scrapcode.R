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
















#=====================================================================================================
#---------------------------Run-Analysis-For-Available-Tissues::ADDIS---------------------------------
#=====================================================================================================
#set.seed(566)

source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")

#check all M1's for CellsEBVtransformedlymphocytes
M1.lymphocytes=ADDIS.M1.check("CellsEBVtransformedlymphocytes")

indexes=c(M1.lymphocytes$Catalog$tiss$CellsEBVtransformedlymphocytes$ADDIS$type2)

triosM1.hic.result=interaction_check(hic.filename="/mnt/ceph/jarredk/HiC_Analyses/Data/lymphoblastoid_cells/ENCFF355OWW.hic", 
                                     trios=indexes[1:10],
                                     resolution=10000,
                                     search.size=100000,
                                     tiss="CellsEBVtransformedlymphocytes",
                                     verbose=TRUE,
                                     plot.h=FALSE,
                                     FDR="ADDIS")






trio.attr=get_trio_attr(trio.index=indexes[1], tissue.name = "CellsEBVtransformedlymphocytes")
cis.data=trio.attr$Attributes$cis
trans.data=trio.attr$Attributes$trans
i=6
search.size=100000
resolution=10000
hic.filename="/mnt/ceph/jarredk/HiC_Analyses/Data/lymphoblastoid_cells/ENCFF355OWW.hic"
#establish snp - gene search "box"
left.pos=paste(cis.data$chr[i],":",
               cis.data$variant_pos[i]-search.size,":",
               cis.data$variant_pos[i]+search.size, sep = "")

#establish trans gene search "box"
right.pos=paste(trans.data$chr[i],":",
                trans.data$left[i]-search.size,":",
                trans.data$right[i]+search.size, sep = "")


data.hic=extract_hic(fileName = hic.filename, 
                     chrs = c(left.pos, right.pos),
                     resol=resolution)

RS=Resample_interactions(filePath = "/mnt/ceph/jarredk/HiC_Analyses/Data/lymphoblastoid_cells/ENCFF355OWW.hic",
                         chrs=c(paste(cis.data$chr[i]),paste(trans.data$chr[i])),
                         res=10000,
                         search.size=search.size,
                         tiss = tiss,
                         trio = trios[i],
                         plot.hist = plot.h,
                         FDR = FDR)
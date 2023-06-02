
library('MRGN', lib = "/mnt/ceph/jarredk/Rpackages/")
tiss=read.csv("/mnt/ceph/jarredk/AddisReRunFiles/tissuenames.csv")
#read in the original master table before changing to protein coding and and lncRNA trios
MT.addis =  loadRData('/mnt/ceph/jarredk/Manuscript/Tables_extra/All_Trios_Master_table_ADDIS.RData')
MT.lond =  loadRData('/mnt/ceph/jarredk/Manuscript/Tables_extra/All_Trios_Master_table_ADDIS.RData')
MT.MRGN = loadRData("/mnt/ceph/jarredk/MRGN_extra/GTEx_tissues_analysis/Protein_LincRNA_PC_Selection/TrioTables/all_master_trio_tables.RData")
bmart=read.table(file="/mnt/ceph/jarredk/Reg_Net/mart_export_merged_lncRNA_fixed.txt", sep = "\t", header = T)
bmart.w.descript = read.csv(file = '/mnt/ceph/jarredk/Manuscript/Biomart_export_with_gene_descriptions.txt')


MT.addis.pr.lnc.only = lapply(c(1:48), function(w,x,y,z) subset(y, Tissue == x[w])[z[[w]]$Trio.orig.index,],
                              x = tiss[,2], y = MT.addis, z = MT.MRGN)


MT.addis.pr.lnc.only2 = lapply(c(1:48), function(x,y,z) y[[x]][,14] = z[[x]]$MRPC.Addis.Inferred.Model,
                               y = MT.addis.pr.lnc.only, z = MT.MRGN)


lapply(MT.addis.pr.lnc.only2, dim)



MT.addis.onetab = do.call('rbind.data.frame', MT.addis.pr.lnc.only2)



MT.onetab = do.call('rbind.data.frame', MT.MRGN)





MRGN.WB.master = loadRData(file = "/mnt/ceph/jarredk/MRGN_extra/GTEx_tissues_analysis/Manuscript_results/trio_tables/all_master_trio_tables.RData")

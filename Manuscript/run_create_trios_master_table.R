

source("/mnt/ceph/jarredk/Manuscript/create_trios_master_table.R")

ms1=create(tissue.names=tissue.names[,2], which.method="ADDIS")
ms2=create(tissue.names=tissue.names[,2], which.method="LOND")

save(ms1, file = "/mnt/ceph/jarredk/Manuscript/All_Trios_Master_table_ADDIS.RData")
save(ms2, file = "/mnt/ceph/jarredk/Manuscript/All_Trios_Master_table_LOND.RData")
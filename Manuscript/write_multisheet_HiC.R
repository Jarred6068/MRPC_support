





tiss=c("CellsCulturedfibroblasts","SkinNotSunExposed","Lung","CellsEBVtransformedlymphocytes")
tisspath=c("fibroblast_cells","Skin","Lung","lymphoblastoid_cells")
tisspath2=c("Fibroblasts","Skin","Lung","lymphocyte")
library(writexl, lib="/mnt/ceph/jarredk/Rpackages")

ADDIS=vector("list", 4)
LOND=vector("list", 4)


names(ADDIS)=tiss
names(LOND)=tiss


for(i in 1:length(tiss)){
  
  ADDIS[[i]]=as.data.frame(read.csv(file = paste0("/mnt/ceph/jarredk/HiC_Analyses/ADDIS/",tisspath[i],"/HiC_",
                                                  tisspath2[i], "_result_ADDIS.csv")))
  
  LOND[[i]]=as.data.frame(read.csv(file = paste0("/mnt/ceph/jarredk/HiC_Analyses/LOND/",tisspath[i],"/HiC_",
                                                 tisspath2[i], "_result_LOND.csv")))

  
  
}


write_xlsx(ADDIS, path = "/mnt/ceph/jarredk/Manuscript/ADDIS_HiC_master.xlsx")
write_xlsx(LOND, path = "/mnt/ceph/jarredk/Manuscript/LOND_HiC_master.xlsx")

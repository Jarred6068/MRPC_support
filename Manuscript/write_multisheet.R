

library('writexl', lib="/mnt/ceph/jarredk/Rpackages")

tiss=read.csv("/mnt/ceph/jarredk/AddisReRunFiles/tissuenames.csv")
library(writexl, lib="/mnt/ceph/jarredk/Rpackages")

AM1T1=vector("list", 48)
AM1T2=vector("list", 48)
LM1T1=vector("list", 48)
LM1T2=vector("list", 48)

names(AM1T1)=tiss[,2]
names(AM1T2)=tiss[,2]
names(LM1T1)=tiss[,2]
names(LM1T2)=tiss[,2]

for(i in 1:dim(tiss)[1]){
  
  AM1T1[[i]]=as.data.frame(read.csv(file = paste0("/mnt/ceph/jarredk/Reg_Net/AM1T1/",tiss[i,2],".csv")))[,-1]
  AM1T2[[i]]=as.data.frame(read.csv(file = paste0("/mnt/ceph/jarredk/Reg_Net/AM1T2/",tiss[i,2],".csv")))[,-1]
  LM1T1[[i]]=as.data.frame(read.csv(file = paste0("/mnt/ceph/jarredk/Reg_Net/LM1T1/",tiss[i,2],".csv")))[,-1]
  LM1T2[[i]]=as.data.frame(read.csv(file = paste0("/mnt/ceph/jarredk/Reg_Net/LM1T2/",tiss[i,2],".csv")))[,-1]
  
  
}


write_xlsx(AM1T1, path = "/mnt/ceph/jarredk/Manuscript/AM1T1_master.xlsx")
write_xlsx(AM1T2, path = "/mnt/ceph/jarredk/Manuscript/AM1T2_master.xlsx")
write_xlsx(LM1T1, path = "/mnt/ceph/jarredk/Manuscript/LM1T1_master.xlsx")
write_xlsx(LM1T2, path = "/mnt/ceph/jarredk/Manuscript/LM1T2_master.xlsx")







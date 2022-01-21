

samples.meta=read.csv("/mnt/ceph/jarredk/Cancer_Selection/Tissue_mut_sampsThu Jan 20 22_02_34 2022.csv")

all.tables=list()
for(i in 1:length(samples.meta$COS.ID.Num)){
  
  all.tables[[i]]=read.csv(paste0("/mnt/ceph/jarredk/Cancer_Selection/genotypes/", 
                                  paste(samples.meta$COS.ID.Num[i], "complexGenotypes.csv.gz",sep="_")))
  
}























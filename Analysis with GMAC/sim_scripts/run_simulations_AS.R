





#---------------------simulation-tables:-------------------------

source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")


#nn=50

#small truth model large inferred model




l1=cross.analyze(tissues="AdiposeSubcutaneous", save=FALSE)

trios=l1$final.tables[[1]]$Trio.Num
med.types=l1$final.tables[[1]]$Mediation.type
print("running Simulations...")

ot1=run.simu12(tissue = "AdiposeSubcutaneous" ,trios=trios, 
               l1.table=l1$final.tables[[1]],
               mod.type.vec=med.types,
               alpha=0.001, n="random")
print("...done")
ot1$Tissue=rep("AdiposeSubcutaneous")


write.table(ot1, file = paste0("/mnt/ceph/jarredk/GMACanalysis/master_tables/simulations_","AdiposeSubcutaneous",".txt"), col.names = T, row.names = F, sep = "\t", 
            quote=F)







#--------------------------------------------------------------------------------------------------

#checking trios with conflicting mediation and permutation p values
#must be run post run_simulations.slurm on cluster
#source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
sim=read.table(file=paste0("/mnt/ceph/jarredk/GMACanalysis/master_tables/simulations_","AdiposeSubcutaneous",".txt"), sep="\t", header = T)

#wt=sim[which(sim$Med.pvalue.GMAC>0.05),]

snp.dist=as.data.frame(matrix(0, nrow = dim(sim)[1], ncol = 3))
colnames(snp.dist)=c("hr","het","ha")

for(i in 1:dim(sim)[1]){
  
  if(sim$Med.type[i]=="Trans.Med"){
    list.data=cross.regress(tissue=sim$Tissue[i], trio.ind=sim$Trio.Num[i], mod.type="trans", addis.pcs=NULL, verbose=F)
  }else{
    list.data=cross.regress(tissue=sim$Tissue[i], trio.ind=sim$Trio.Num[i], mod.type="cis", addis.pcs=NULL, verbose=F)
  }
  
  
  #account for SNPs with less than 3 genotypes
  if(length(summary(factor(list.data$GMAC$SNP)))<3){
    
    if( length(which(names(summary(factor(list.data$GMAC$SNP)))=="0"))==0 ){
      snp.dist[i,]=c(NA ,summary(factor(list.data$GMAC$SNP)))
    }else if ( length(which(names(summary(factor(list.data$GMAC$SNP)))=="1"))==0 ){
      snp.dist[i,]=c(summary(factor(list.data$GMAC$SNP))[1], NA , summary(factor(list.data$GMAC$SNP))[2])
    }else{
      snp.dist[i,]=c(summary(factor(list.data$GMAC$SNP)), NA)
    }
    
  }else{
    
    snp.dist[i,]=summary(factor(list.data$GMAC$SNP))
    
  }
  
}


output1=cbind.data.frame(sim, snp.dist)

write.csv(output1, file = paste0("/mnt/ceph/jarredk/GMACanalysis/master_tables/TRIOS_imbalanced_genotypes_","AS",".csv"), 
          row.names = F)





print("...Finished")



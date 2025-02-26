





#---------------------simulation-tables:-------------------------

source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")


#nn=50

#small truth model large inferred model
tiss=c("WholeBlood","AdiposeSubcutaneous")
tiss.save=c("WB","AS")
for(i in 1:2){



   l1=cross.analyze(tissues=tiss[i], save=FALSE)

   trios=l1$final.tables[[1]]$Trio.Num
   med.types=l1$final.tables[[1]]$Mediation.type
  print("running Simulations...")

  ot1=run.simu12(tissue = tiss[i] ,trios=trios, 
                 l1.table=l1$final.tables[[1]],
                 mod.type.vec=med.types,
                 alpha=0.001, n="random")
  print("...done")
  ot1$Tissue=rep(tiss[i])
  

    write.table(ot1, file = paste0("/mnt/ceph/jarredk/GMACanalysis/master_tables/simulations_",tiss[i],".txt"), col.names = T, row.names = F, sep = "\t", 
                quote=F)



    
    
    #--------------------------------------------------------------------------------------------------
    
    #Add in simulation for MRPC on GMAC model to the simulations.txt dataset
    #source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
    sim=read.table(file=paste0("/mnt/ceph/jarredk/GMACanalysis/master_tables/simulations_",tiss[i],".txt"), sep="\t", header = T)
    
    p=NULL
    
    for(i in 1:dim(sim)[1]){
      print(i)
      
      L2A=Lond2Addis.lookup(trio.index=sim$Trio.Num[i], tissue.name=sim$Tissue[i], with.pc=TRUE)
      
      if(length(colnames(L2A$correlation))>3){
        addis.pcs=colnames(L2A$correlation)[-c(1:3)]
      }else{
        addis.pcs=NULL
      }
      
      list.data=cross.regress(tissue=sim$Tissue[i], trio.ind=sim$Trio.Num[i], mod.type=sim$Med.type[i], addis.pcs=addis.pcs, verbose=F)
      
      p[i]=simu3(MRPC.data=list.data$addis, GMAC.data=list.data$GMAC, mod.type=sim$Med.type[i], verbose=F)$pvalue
      
      
    }
    
    sim$TGM.p=p
    
    write.table(sim, file = paste0("/mnt/ceph/jarredk/GMACanalysis/master_tables/simulations_",tiss[i],".txt"), 
                col.names = T, row.names = F, sep = "\t", quote=F)
    
    
    #--------------------------------------------------------------------------------------------------
    
    #checking trios with conflicting mediation and permutation p values
    #must be run post run_simulations.slurm on cluster
    #source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
    sim=read.table(file=paste0("/mnt/ceph/jarredk/GMACanalysis/master_tables/simulations_",tiss[i],".txt"), sep="\t", header = T)
    
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
    
    write.csv(output1, file = paste0("/mnt/ceph/jarredk/GMACanalysis/master_tables/TRIOS_imbalanced_genotypes_",tiss.save[i],".csv"), 
              row.names = F)
    
    
    #-------------------------------------------------------------------------------------------------
    
    #source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
    #Add in permutation pvalue for simulation for MRPC on GMAC model to the simulations.txt dataset
    t1=read.csv(file = paste0("/mnt/ceph/jarredk/GMACanalysis/master_tables/TRIOS_imbalanced_genotypes_",tiss.save[i],".csv"))
    
    p=NULL
    
    
    for(i in 1:dim(t1)[1]){
      #print(i)
      
      L2A=Lond2Addis.lookup(trio.index=t1$Trio.Num[i], tissue.name=t1$Tissue[i], with.pc=TRUE)
      print(dim(L2A$correlation))
      if(length(colnames(L2A$correlation))>3){
        addis.pcs=colnames(L2A$correlation)[-c(1:3)]
      }else{
        addis.pcs=NULL
      }
      
      list.data=cross.regress(tissue=t1$Tissue[i], trio.ind=t1$Trio.Num[i], mod.type=t1$Med.type[i], addis.pcs=addis.pcs, verbose=F)
      
      out=simu3(MRPC.data=list.data$addis, GMAC.data=list.data$GMAC, mod.type=t1$Med.type[i], verbose=F)
      print(dim(out$sim.data))
      p[i]=run.permuted.reg(out$sim.data, nperms=1000, Alg="ADDIS", med.type = t1$Med.type[i])$p.value
      print(p[i])
      
    }
    
    
    t1$TGM.perm=p
    
    
    write.csv(t1, file = paste0("/mnt/ceph/jarredk/GMACanalysis/master_tables/TRIOS_imbalanced_genotypes_",tiss.save[i],".csv"), 
              row.names = F)
    
}

print("...Finished")



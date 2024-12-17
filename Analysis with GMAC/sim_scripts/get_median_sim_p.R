#code to run multiple simulations to check the stability of the analysis

#load in needed stuff
sim=read.csv(file="/mnt/ceph/jarredk/GMACanalysis/master_tables/TRIOS_imbalanced_genotypes_WB.csv", header = T)
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")

tissue="WholeBlood"
l1=cross.analyze(tissues=tissue, save=FALSE)
#pull out trios which were no longer significant after the GSS simulation 
#trios=sim$Trio.Num[sim$Nominal.p.GSS>0.05][1:5]
#mod.type.vec=sim$Med.type[sim$Nominal.p.GSS>0.05][1:5]
trios = sim$Trio.Num
mod.type.vec = sim$Med.type

ss=1000
storage1=vector("list", length = length(trios))
plot.it=FALSE
STM.median.p=NULL
LTM.median.p=NULL
TGM.median.p=NULL
GSS.median.p=NULL

for(i in 1:length(trios)){
  
  
  print(paste0("getting data for trio # ", trios[i]))
  if(mod.type.vec[i]=="Both"){
    
    list.data=cross.regress(tissue=tissue, 
                            trio.ind=trios[i], 
                            mod.type="cis", 
                            #addis.pcs=addis.pcs, 
                            verbose = FALSE)
    
  }else if(mod.type.vec[i]=="Cis.Med"){
    
    list.data=cross.regress(tissue=tissue, 
                            trio.ind=trios[i], 
                            mod.type="cis", 
                            #addis.pcs=addis.pcs, 
                            verbose = FALSE)
  }else{
    
    list.data=cross.regress(tissue=tissue, 
                            trio.ind=trios[i], 
                            mod.type="trans", 
                            #addis.pcs=addis.pcs, 
                            verbose = FALSE)
    
  }
  print("...done...")
  
  p1=as.data.frame(matrix(0, nrow = ss, ncol = 4))
  b1=as.data.frame(matrix(0, nrow = ss, ncol = 4))
  colnames(p1)=c("simu1.p", "simu2.p", "simu3.p", "simu4.p")
  colnames(b1)=c("simu1.b", "simu2.b", "simu3.b", "simu4.b")
  
  print(paste0("running all simulations for trio # ", trios[i]))
  for(j in 1:ss){
    
    
    out=simu1(list.data$GMAC, alpha = 0.001, mod.type=mod.type.vec[i], verbose=F)
    nn=floor(runif(1, 1, 20))
    #print(nn)
    out2=simu2(tissue = tissue, data=list.data$GMAC, mod.type=mod.type.vec[i], n=nn, verbose=F)
    out3=simu3(MRPC.data=list.data$addis, GMAC.data=list.data$GMAC, mod.type=mod.type.vec[i], verbose=F)
    out4=simu4(GMAC.data=list.data$GMAC, mod.type=mod.type.vec[i], verbose=F)
    
    p1$simu1.p[j] = out$pvalue
    p1$simu2.p[j] = out2$pvalue
    p1$simu3.p[j] = out3$pvalue
    p1$simu4.p[j] = out4$pvalue
    
    b1$simu1.b[j] = out$b.coef
    b1$simu2.b[j] = out2$b.coef
    b1$simu3.b[j] = out3$b.coef
    b1$simu4.b[j] = out4$b.coef
    print(p1$simu1.p)
  }
  print("...done!")
  print(head(p1))
  
  STM.median.p[i] = median(p1$simu1.p)
  LTM.median.p[i] = median(p1$simu2.p)
  TGM.median.p[i] = median(p1$simu3.p)
  GSS.median.p[i] = median(p1$simu4.p)
  
  
  storage1[[i]] = cbind.data.frame(p1, b1)
  
  
  if(plot.it==TRUE){
    png(paste0("/mnt/ceph/jarredk/GMACanalysis/Unstable_trios_plots/plot_", tissue, "_", trios[i], "_",mod.type.vec[i], ".png"))
    plot(p1, p3, pch=21, bg="black",
         xlab="Simulated Para. P-value",
         ylab="Simulated Nominal P-value",
         main = paste0("GSS 1000 simulations of trio = ", trios[i], ";", tissue ))
    abline(a=0, b=1, lty="dotted", col="red")
    dev.off()
    
    png(paste0("/mnt/ceph/jarredk/GMACanalysis/Unstable_trios_plots/histogram_", tissue, "_", trios[i], "_",mod.type.vec[i], ".png"))
    par(mfrow=c(2,2))
    hist(p1, xlab="Simulated Para. P-value",
         ylab="Simulated Nominal P-value",
         main = paste0("GSS 1000 simulations of trio = ", trios[i], ";", tissue ))
    hist(p2,
         xlab="Simulated Para. P-value",
         ylab="Simulated Nominal P-value",
         main = paste0("GSS 1000 simulations of trio = ", trios[i], ";", tissue ))
    hist(p3,
         xlab="Simulated Para. P-value",
         ylab="Simulated Nominal P-value",
         main = paste0("GSS 1000 simulations of trio = ", trios[i], ";", tissue ))
    dev.off()
  }
}


sim.new=cbind.data.frame(sim, STM.median.p, LTM.median.p, TGM.median.p, GSS.median.p)

save(storage1, file = "/mnt/ceph/jarredk/GMACanalysis/master_tables/multi_sim_data.RData")
write.csv(sim.new, file = "/mnt/ceph/jarredk/GMACanalysis/master_tables/TRIOS_imbalanced_genotypes_WB_Updated_12_10_2021.csv")


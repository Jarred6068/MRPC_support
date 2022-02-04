
#load in a few useful functions from previous coding 
source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")
#load M0
m0.regres=loadRData(fileName = "/mnt/ceph/jarredk/MRPC_UPDATE/Simulations/Sim_result_Data/M0.reg.res.RData")
m0.inf=loadRData(fileName = "/mnt/ceph/jarredk/MRPC_UPDATE/Simulations/Sim_result_Data/M0.inf.mods.RData")
acc.m0=sapply(m0.inf, FUN = function(x) length(which(x=="M0.1" | x=="M0.2"))/length(x) )
#load M1
m1.regres=loadRData(fileName = "/mnt/ceph/jarredk/MRPC_UPDATE/Simulations/Sim_result_Data/M1.reg.res.RData")
m1.inf=loadRData(fileName = "/mnt/ceph/jarredk/MRPC_UPDATE/Simulations/Sim_result_Data/M1.inf.mods.RData")
acc.m1=sapply(m1.inf, FUN = function(x) length(which(x=="M1.1" | x=="M1.2"))/length(x) )
#load M2
m2.regres=loadRData(fileName = "/mnt/ceph/jarredk/MRPC_UPDATE/Simulations/Sim_result_Data/M2.reg.res.RData")
m2.inf=loadRData(fileName = "/mnt/ceph/jarredk/MRPC_UPDATE/Simulations/Sim_result_Data/M2.inf.mods.RData")
acc.m2=sapply(m2.inf, FUN = function(x) length(which(x=="M2.1" | x=="M2.2"))/length(x) )
#load M3
m3.regres=loadRData(fileName = "/mnt/ceph/jarredk/MRPC_UPDATE/Simulations/Sim_result_Data/M3.reg.res.RData")
m3.inf=loadRData(fileName = "/mnt/ceph/jarredk/MRPC_UPDATE/Simulations/Sim_result_Data/M3.inf.mods.RData")
acc.m3=sapply(m3.inf, FUN = function(x) length(which(x=="M3"))/length(x) )
#load M4
m4.regres=loadRData(fileName = "/mnt/ceph/jarredk/MRPC_UPDATE/Simulations/Sim_result_Data/M4.reg.res.RData")
m4.inf=loadRData(fileName = "/mnt/ceph/jarredk/MRPC_UPDATE/Simulations/Sim_result_Data/M4.inf.mods.RData")
acc.m4=sapply(m4.inf, FUN = function(x) length(which(x=="M4"))/length(x) )

mod.inf.list=list(m0.inf,m1.inf,m2.inf,m3.inf,m4.inf)
acc.list=list(acc.m0,acc.m1,acc.m2,acc.m3,acc.m4)
save(acc.list, file = "/mnt/ceph/jarredk/MRPC_UPDATE/Simulations/Sim_result_Data/all.acc.noRARE.RData")


#get the accuracy list from local machine
load("/mnt/ceph/jarredk/MRPC_UPDATE/Simulations/Sim_result_Data/all.acc.noRARE.RData")
#load("C:/Users/Bruin/Documents/GitHub/MRPC_support/MRPC_UPDATE/Simulations/all.acc.noRARE.RData")
#set initial parameters for simulation
#sample sizes
n=c(50, 100, 500, 1000)
#noise in the data/SD of errors
noise=c(0.2, 0.5, 0.8, 1.2)
#frequency of the minor allele
minor.allele=c(0.1, 0.2, 0.3, 0.4)
#signal strength of edge
b1.1=c(0.2, 0.4, 0.6, 0.8)
b1.2=c(0.2, 0.4, 0.6, 0.8)

#expand to all get all combos of simulation conditions
model.params=expand.grid(n, noise, minor.allele, b1.1, b1.2)
colnames(model.params)=c("sample.size","SD", "minor.freq", "b1.1","b1.2")

get.all.acc=function(acc=acc.list, params.map=model.params, plot=TRUE){

  plot.acc=vector("list",length=5)
  names(plot.acc)=c("M0","M1","M2","M3","M4")
  
  for(i in 1:5){
    
    mean.acc=as.data.frame(matrix(0, nrow = 4, ncol = 5))
    colnames(mean.acc)=c("sample.size","SD", "minor.freq", "b1.1","b1.2")
    sim.param=as.data.frame(matrix(0, nrow = 4, ncol = 5))
    colnames(sim.param)=c("sample.size","SD", "minor.freq", "b1.1","b1.2")
    
    for(j in 1:5){
      
      for(k in 1:4){
        #average accruacy across sample size
        mean.acc[k,j]=mean(acc[[i]][which(params.map[,j]==unique(params.map[,j])[k])])
        sim.param[k,j]=unique(params.map[,j])[k]
      }
    }
    
    plot.acc[[i]]=mean.acc
    #print(sim.param)
  }
  
  if(plot==T){
    xlabels=c("sample.size","SD", "minor.freq", "b1.1","b1.2")
    ylabel="Mean Accuracy"
    main.labels=paste0("Mean Acc Across ",xlabels)
    colors1=c("black","red","blue","green2","yellow3")
    xlim.upper=c(1100, 1.4, 0.5, 1, 1)
    
    pdf("/mnt/ceph/jarredk/MRPC_UPDATE/Simulations/all.mods.plots.pdf")
    
    par(mfrow=c(3,2), xpd=TRUE)
    for(p in 1:5){
      
      plot(x=sim.param[,p], 
           y=plot.acc[[1]][,p], 
           type = "b", 
           pch=21, 
           col=colors1[1], 
           bg=colors1[1],
           lty=1,
           lwd=2,
           bty="L",
           main = main.labels[p],
           ylab = ylabel,
           xlab = xlabels[p],
           xlim = c(min(sim.param[,p]),xlim.upper[p]),
           ylim = c(0,1))
      # legend("bottomleft", 
      #        inset=c(-0.2,0), 
      #        col=colors1, 
      #        legend = c("M0","M1","M2","M3","M4"), 
      #        pch = rep(21,5), 
      #        fill = colors1,
      #        horiz = T)
      for(i in 2:5){
        
        lines(x=sim.param[,p], 
              y=plot.acc[[i]][,p],
              type = "b",
              pch=21,
              lty=i,
              lwd=2,
              col=colors1[i],
              bg=colors1[i])
      }
      text(x=rep(xlim.upper[p],5), 
           y=seq(0.3, 0.9, 0.15), 
           labels = c("M0","M1","M2","M3","M4"),
           col = colors1)
    }
    par(mfrow=c(1,1))
    dev.off()
  }
  
  
  return(list(acc=plot.acc, params=sim.param))
  
}

mean.acc=get.all.acc()























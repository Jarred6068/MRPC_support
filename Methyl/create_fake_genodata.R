
set.seed(123)

fake_genos=function(csize=NULL, rsize=NULL){
  
  genos=sample(c(0,1,2), csize*rsize, replace = TRUE)
  
  genos.mat=as.data.frame(matrix(genos, nrow = rsize, ncol = csize, byrow = T))
  
  return(genos.mat)
  
}

genos.mat=fake_genos(4000000, 1000)

genos.mat[1:10,1:5]

sim.chr=sample(c(1:22,"X","Y"), dim(genos.mat)[2], replace = TRUE)

simu.pos=function(genosmat=NULL, chr.vec=NULL){
  
  #pos=matrix(0, nrow = dim(genosmat)[2], ncol = 2)
  pos=rep(0, 4000000)
  chr.sizes=read.csv(file = "/mnt/ceph/jarredk/Methyl/chrsizes.csv")
  
  for(i in 1:length(chr.vec)){
    
    idx=match(chr.vec[i],chr.sizes[,1])
    start.pos=runif(1, 1000000, chr.sizes[idx,2]-1000000)
    #end.pos=start.pos+1000000
    #print(start.pos)
    #print(end.pos)
    pos[i]=start.pos
    
  }
  
  return(pos)
  
}

names1=paste0("SNP","_",c(1:4000000),"_", sim.chr)

sim.chrpos=simu.pos(genosmat = genos.mat, chr.vec = sim.chr)
simulated.meta=cbind.data.frame(sim.chr, sim.chrpos)

colnames(genos.mat)=names1

genos.mat[1:10,1:5]
sim.chr[1:5]

save(simulated.meta, file = "/mnt/ceph/jarredk/Methyl/fakegenoMeta_2.Rdata")
#simulated.meta=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/fakegenoMeta.Rdata")

save(genos.mat, file = "/mnt/ceph/jarredk/Methyl/fakegenosBIG_2.Rdata")

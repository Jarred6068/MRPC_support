
set.seed(123)

n=100
m=1000

fake_genos=function(csize=NULL, rsize=NULL){
  
  genos=sample(c(0,1,2), csize*rsize, replace = TRUE)
  
  genos.mat=as.data.frame(matrix(genos, nrow = rsize, ncol = csize, byrow = T))
  
  return(genos.mat)
  
}

genos.mat=fake_genos(n, m)

genos.mat[1:10,1:5]

sim.chr=sample(c(1:22,"X","Y"), dim(genos.mat)[2], replace = TRUE)

simu.pos=function(genosmat=NULL, chr.vec=NULL, q=n){
  
  #pos=matrix(0, nrow = dim(genosmat)[2], ncol = 2)
  pos=rep(0, q)
  chr.sizes=read.csv(file = "/mnt/ceph/jarredk/Methyl/chrsizes.csv")
  
  for(i in 1:length(chr.vec)){
    
    idx=match(chr.vec[i],chr.sizes[,1])
    start.pos=round(runif(1, 1000000, chr.sizes[idx,2]-1000000))
    #end.pos=start.pos+1000000
    #print(start.pos)
    #print(end.pos)
    pos[i]=start.pos
    
  }
  
  return(pos)
  
}
possible.alleles=c("A","T","C","G")

sim.other=sample(possible.alleles, dim(genos.mat)[2], replace = TRUE)
sim.counted=NULL

fn=function(x,v){x=sample(v[-which(x==v)], 1, replace = T )}

sim.counted=unlist(lapply(as.list(sim.other), fn, possible.alleles))

sim.chrpos=simu.pos(genosmat = genos.mat, chr.vec = sim.chr)

names1=paste0("SNP","_",sim.chrpos,"_chr", sim.chr)
simulated.meta=cbind.data.frame(names1, sim.chr, sim.chrpos, sim.other, sim.counted)

#colnames(genos.mat)=names1
colnames(simulated.meta)=c("SNP_ID","chr","coordinate", "OtherAllele", "CountedAllele")

genos.mat=t(genos.mat)
genos.mat.new=cbind.data.frame(names1, genos.mat)
colnames(genos.mat.new)=c("SNP_ID", paste0("bsgs_", c(1:1000)))

genos.mat.new[1:10,1:5]
sim.chr[1:5]

#save(simulated.meta, file = "/mnt/ceph/jarredk/Methyl/fakegenoMeta_2.Rdata")
write.table(simulated.meta, file = "/mnt/ceph/jarredk/Methyl/CC_all/Correlation_CalculationV5_Scott/Input/Example_Genotype_Metadata_small.txt", 
            col.names = T, sep = "\t", quote = F, row.names = F)

#simulated.meta=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/fakegenoMeta.Rdata")

#save(genos.mat.new, file = "/mnt/ceph/jarredk/Methyl/fakegenosBIG_2.Rdata")

write.table(genos.mat.new, file = "/mnt/ceph/jarredk/Methyl/CC_all/Correlation_CalculationV5_Scott/Input/Example_Genotype_data_small.txt", 
            col.names = T, sep = "\t", quote = F, row.names = F)

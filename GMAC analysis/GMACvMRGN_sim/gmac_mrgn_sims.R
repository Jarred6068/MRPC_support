
library("MRGN", lib="/mnt/ceph/jarredk/Rpackages")
library("GMAC", lib="/mnt/ceph/jarredk/Rpackages")
tissue.name="WholeBlood"

sim.data=loadRData("/mnt/ceph/jarredk/MRGN_extra/Simulations/Sim_wPC_result_data/mrgn_10k_datasets_wpc2small_all_mods.RData")
num.trios=length(sim.data)
sim.data.red=lapply(sim.data, function(x)x[,1:3])
sim.data.mat=do.call("cbind", sim.data.red)
snp.idx=seq(1, dim(sim.data.mat)[2]-2, 3)
snp.dat.cis=sim.data.mat[,snp.idx]
exp.dat=sim.data.mat[,-snp.idx]
trios.idx=cbind(c(1:num.trios), matrix(c(1:dim(exp.dat)[2]), nrow = num.trios, ncol = 2, byrow = T))
pc.matrix=loadRData(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissue.name,
                          "_AllPC/PCs.matrix.",tissue.name,".RData", sep = ""))



dim(t(exp.dat))
dim(t(snp.dat.cis))
dim(t(pc.matrix))

input.list=list(cov.pool=t(pc.matrix), snp.dat.cis=t(snp.dat.cis), exp.dat=t(exp.dat), trio.indexes=trios.idx,
                nperm=10000, clusters=cl, use.nominal.p=TRUE)

cl <- makeCluster(10)


output.cis.med <- gmac(cl=input.list$clusters, known.conf = input.list$known.conf, cov.pool = input.list$cov.pool,
                       exp.dat = input.list$exp.dat, snp.dat.cis = input.list$snp.dat.cis,
                       trios.idx = input.list$trio.indexes, nperm = input.list$nperm,
                       nominal.p = input.list$use.nominal.p)


stopCluster(cl)


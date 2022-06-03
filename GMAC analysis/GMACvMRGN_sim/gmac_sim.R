
library("GMAC", lib="/mnt/ceph/jarredk/Rpackages")
#cl <- makeCluster(10)
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")
library("MRGN",lib="/mnt/ceph/jarredk/Rpackages")
library("qvalue", lib="/mnt/ceph/jarredk/Rpackages")

small.datasets=loadRData(file = "/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/simulated_data/mrgn_v_gmac_v_mrpc_10k_datasets_c_kc_all_mods.RData")
kc = t(as.matrix(output.WB$input.list$known.conf))
#load in pc mat from whole blood
tissue.name="WholeBlood"
pc.matrix=loadRData(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissue.name,
                          "_AllPC/PCs.matrix.",tissue.name,".RData", sep = ""))
#perform GMAC inferences
#GMAC
print("Running GMAC.....")
num.trios=length(small.datasets)
small.datasets.red=lapply(small.datasets, function(x)x[,1:3])
small.datasets.mat=do.call("cbind", small.datasets.red)
snp.idx=seq(1, dim(small.datasets.mat)[2]-2, 3)
snp.dat.cis=small.datasets.mat[,snp.idx]
exp.dat=small.datasets.mat[,-snp.idx]
trios.idx=cbind(c(1:num.trios), matrix(c(1:dim(exp.dat)[2]), nrow = num.trios, ncol = 2, byrow = T))


dim(t(exp.dat))
dim(t(snp.dat.cis))
dim(t(pc.matrix))

input.list=list(kc = t(kc), cov.pool=t(pc.matrix), snp.dat.cis=t(snp.dat.cis), exp.dat=t(exp.dat), trio.indexes=trios.idx,
                nperm=1000, use.nominal.p=TRUE)



output.cis.med <- gmac(known.conf = input.list$kc, cov.pool = input.list$cov.pool,
                       exp.dat = input.list$exp.dat, snp.dat.cis = input.list$snp.dat.cis,
                       trios.idx = input.list$trio.indexes, nperm = input.list$nperm,
                       nominal.p = input.list$use.nominal.p)

print(output.cis.med)

output.trans.med <- gmac(known.conf = input.list$kc, cov.pool = input.list$cov.pool,
                         exp.dat = input.list$exp.dat, snp.dat.cis = input.list$snp.dat.cis,
                         trios.idx = input.list$trio.indexes[,c(1,3,2)], nperm = input.list$nperm,
                         nominal.p = input.list$use.nominal.p)
#stopCluster(cl)
print(output.trans.med)

#reorganize output for saving
#cis results
out.table.cis=cbind.data.frame(output.cis.med[[1]], output.cis.med[[2]])
colnames(out.table.cis)=c(paste0('pval_', colnames(output.cis.med[[1]])),
                          paste0('effect_change_', colnames(output.cis.med[[2]])))

out.list.cis=list(out.table.cis, input.list, output.cis.med[[3]])
names(out.list.cis)=c("output.table", "input.list", "cov.indicator.list")

#trans results
out.table.trans=cbind.data.frame(output.trans.med[[1]], output.trans.med[[2]])
colnames(out.table.cis)=c(paste0('pval_', colnames(output.trans.med[[1]])),
                          paste0('effect_change_', colnames(output.trans.med[[2]])))

out.list.trans=list(out.table.trans, input.list, output.trans.med[[3]])
names(out.list.trans)=c("output.table", "input.list", "cov.indicator.list")



print("Done!...Saving...")

save(out.list.trans, file = "/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/simulated_data/gmac_10k_trans_results_c_kc_all_mods.RData")
save(out.list.cis, file = "/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/simulated_data/gmac_10k_cis_results_c_kc_all_mods.RData")


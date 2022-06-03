

library("MRGN",lib="/mnt/ceph/jarredk/Rpackages")
source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")

small.datasets=loadRData(file = "/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/simulated_data/mrgn_v_gmac_v_mrpc_10k_datasets_c_kc_all_mods.RData")
just.trios=lapply(small.datasets, function(x) x[,1:3])
tissue.name="WholeBlood"
pc.matrix=loadRData(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissue.name,
                          "_AllPC/PCs.matrix.",tissue.name,".RData", sep = ""))
confs=get.conf(trios=just.trios, PCscores=pc.matrix, blocksize = 2000, method = "correlation")

conf.list=lapply(confs$sig.asso.pcs, function(x,y){y[,x[[1]]]}, y=pc.matrix)

trios.with.pcs=mapply(cbind.data.frame, just.trios, conf.list)

#add in known confounders

kc = t(as.matrix(output.WB$input.list$known.conf))

trios.with.pcs2=lapply(trios.with.pcs, function(x,y) cbind.data.frame(x,y), y = kc)

print(lapply(trios.with.pcs2[1:5], head))

#preform regressions and classify model types
#MRGN
print("Running MRGN")
reg.res=NULL
inf.mods=NULL
#regression
reg.res=sapply(trios.with.pcs2, infer.trio, nperm=1000)
#model class
inf.mods=apply(reg.res, 2, class.vec)

#get estimate of time to compute each trio
# for(i in 1:length(small.datasets)){
#   start.time=Sys.time()
#   z=infer.trio(small.datasets[[i]])
#   x=class.vec(z)
#   end.time=Sys.time()
#   params$Time.to.compute.mrgn[i]=difftime(end.time, start.time, units = 'mins')
# }
print("Done!...Saving results...")

#save
save(reg.res, file = "/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/simulated_data/mrgn_10k_regres_results_c_kc_all_mods.RData")
save(inf.mods, file = "/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/simulated_data/mrgn_10k_inf_results_c_kc_all_mods.RData")

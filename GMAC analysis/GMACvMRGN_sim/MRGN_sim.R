

library("MRGN",lib="/mnt/ceph/jarredk/Rpackages")

small.datasets=loadRData(file = "/mnt/ceph/jarredk/GMACanalysis/GMACvMRGN_sim/simulated_data/mrgn_v_gmac_v_mrpc_10k_datasets_c_kc_all_mods.RData")
#preform regressions and classify model types
#MRGN
print("Running MRGN")
reg.res=NULL
inf.mods=NULL
#regression
reg.res=sapply(small.datasets, infer.trio, nperm=1000)
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

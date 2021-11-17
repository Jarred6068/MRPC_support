#script to run GMACanalysisV2

library('GMAC', lib='/mnt/ceph/jarredk/Rpackages')
library('mice', lib='/mnt/ceph/jarredk/Rpackages')
library('missMDA', lib='/mnt/ceph/jarredk/Rpackages')

source("/mnt/ceph/jarredk/GMACanalysis/GMACanalysis.R")

run.GMAC(tissues.vec=tissues.vec[2, ], path.tables=path, mediation.type='cis')

run.GMAC(tissues.vec=tissues.vec[2, ], path.tables=path, mediation.type='trans')
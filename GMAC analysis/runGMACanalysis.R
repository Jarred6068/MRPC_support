
library('GMAC', lib='/mnt/ceph/jarredk/Rpackages')
library('mice', lib='/mnt/ceph/jarredk/Rpackages')
library('missMDA', lib='/mnt/ceph/jarredk/Rpackages')

source("/mnt/ceph/jarredk/GMACanalysis/GMACanalysisV2.R")

run.GMAC(tissues.vec=tissues.vec, path.tables=path)


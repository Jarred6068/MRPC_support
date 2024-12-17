#script to run GMACanalysis

library('GMAC', lib='/mnt/ceph/jarredk/Rpackages')
library('mice', lib='/mnt/ceph/jarredk/Rpackages')
library('missMDA', lib='/mnt/ceph/jarredk/Rpackages')

source("/mnt/ceph/jarredk/GMACanalysis/GMACanalysis.R")

source("/mnt/ceph/jarredk/GMACanalysis/GMACpostproc.R")

l1=cross.analyze(tissues="AdiposeSubcutaneous", save=FALSE)
l1.tab=l1$final.tables[[1]]

trans.trios=l1.tab$Trio.Num[which(l1.tab$Mediation.type=="Trans.Med")]
cis.trios=l1.tab$Trio.Num[which(l1.tab$Mediation.type!="Trans.Med")]

run.GMAC(tissues.vec=tissues.vec[1,], path.tables=path, mediation.type='cis', which.trios=cis.trios)

run.GMAC(tissues.vec=tissues.vec[1,], path.tables=path, mediation.type='trans', which.trios=trans.trios)
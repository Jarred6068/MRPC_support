

source("/mnt/ceph/jarredk/MRPC_UPDATE/GTEx_tissues_analysis/triospecific_GTEx.R")

result=mrnet.analyze(GTEx.tissname="WholeBlood", save=TRUE)


pvals=PermReg(trio = A, 
              t.obs21 = pt.out$tvals[2], 
              t.obs22 = pt.out$tvals[4], 
              p11 = pt.out$pvals[1],
              p12 = pt.out$pvals[3],
              m = 1000)






























AA.sum=read.csv(file='C:/Users/Bruin/Desktop/Research Assistantship/ADDIS_Analysis_summary.csv', header = TRUE)


idx1=seq(6,16,2)
idx2=seq(7,17,2)
head(AA.sum[,idx1])
head(AA.sum[,idx2])


#calculate the difference in LOND to ADDIS as val of M(i) from ADDIS - M(i) from LOND
diffs=as.data.frame(matrix(0, nrow = dim(AA.sum)[1], ncol = length(idx1) ))
for(i in 1:length(idx1)){
  
  diffs[,i]=AA.sum[,idx1[i]]-AA.sum[,idx2[i]]
  
}

namer=c("M0","M1","M2","M3","M4","Other")

colnames(diffs)=namer
row.names(diffs)=as.character(AA.sum[,1])
print(diffs)

#check to make sure sample sizes match:

AA.sum$X..Trios.tested-rowSums(AA.sum[,idx1])

#================================Check__1:WholeBlood==============================================
#starting with WholeBlood:
#find where M1 models classified under LOND went...
#library(rlist, lib="/mnt/ceph/jarredk/Rpackages")
list.LOND=loadRData(file="/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/WholeBlood_AllPC/List.models.WholeBlood.all.RData")
list.ADDIS=loadRData(file="/mnt/ceph/jarredk/AddisReRunFiles/List.models.WholeBlood.all.RData")

summary(list.LOND)
summary(list.ADDIS)


matched.trios=match(list.LOND$M0, list.ADDIS$M0)
matched=list.LOND$M1[!is.na(matched.trios)]
notmatched=list.LOND$M1[is.na(matched.trios)]



location.ADDIS=list()
for(i in 1:length(list.LOND)){
  
  indicies=match(notmatched, list.ADDIS[[i]])
  indicies=indicies[!is.na(indicies)]
  print(indicies)
  location.ADDIS[[i]]=list.ADDIS[[i]][indicies]
  
}

names(location.ADDIS)=c("match.MO","match.M1","match.M2","match.M3","match.M4")


print(notmatched)

print(location.ADDIS)





#================================Check__1:NerveTibial==============================================


load(file="/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/NerveTibial_AllPC/List.models.NerveTibial.all.RData")
list.LOND=List.models.NerveTibial.all
load(file="/mnt/ceph/jarredk/RtissuesData/NerveTibial/List.models.NerveTibial.all.RData")
list.ADDIS=List.models.NerveTibial.all

summary(list.LOND)
summary(list.ADDIS)


sim.trios=match(list.LOND$M1, list.ADDIS$M1)
sim.trios
nomatch=list.LOND$M1[is.na(sim.trios)]
nomatch
length(nomatch)


location.ADDIS=list()
for(i in 1:length(list.LOND)){
  
  indicies=match(nomatch, list.ADDIS[[i]])
  indicies=indicies[!is.na(indicies)]
  print(indicies)
  location.ADDIS[[i]]=list.ADDIS[[i]][indicies]
  
}

names(location.ADDIS)=c("match.MO","match.M1","match.M2","match.M3","match.M4")


print(nomatch)

print(location.ADDIS)















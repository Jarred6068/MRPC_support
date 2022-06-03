

source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")

#--------------------helper_fn------------------------

find.models=function(list.of.models=NULL, vec.size=NULL){

  models=c("M0","M1","M2","M3","M4")

  chr.vec=rep(NA, vec.size)

  for(i in 1:length(list.of.models)){

    chr.vec[list.of.models[[i]]]=models[i]
  }

  chr.vec[which(is.na(chr.vec))]="Other"

  return(chr.vec)
}

#--------------------helper_fn------------------------

#
# match.fn=function(x, indata){
#
#   return(indata[which(indata$Stable.ID==x),])
#
# }

#--------------------helper_fn------------------------

replace.vals=function(x){

  x=replace(x, which(x==""), NA)

  return(x)
}



#--------------------helper_fn------------------------


fill.table=function(X, cis=NULL, trans=NULL, size=NULL){


  #fill in cis data
  for(i in 1:dim(cis)[1]){

    X[i,c(3,6,8:9,12)]=cis[i,c(4:5,1:2,3)]
    X[i,c(5,7,10:11,13)]=trans[i,c(4:5,1:2,3)]
  }

  return(X)

}


#--------------------MAIN_fn------------------------
create=function(tissue.names=tissue.names[,2], which.method="ADDIS", save.table=TRUE){


  #loading bar 2
  print("---Creating Master Table...---")
  lbar2 <- txtProgressBar(min = 0, max = length(tissue.names), style = 3)

  filen="/mnt/ceph/jarredk/Reg_Net/mart_export_merged_lncRNA_fixed.txt"
  bm.data=read.table(file=filen, sep="\t", header=T)


  bm.data=as.data.frame(apply(bm.data, 2, replace.vals))
  colnames(bm.data)=c("Gene.start", "Gene.end", "Type", "Name", "Chr", "Stable.ID")

  master.list=vector("list", length = length(tissue.names))

  for(t in 1:length(tissue.names)){

    Sys.sleep(0.005)
    setTxtProgressBar(lbar2, t)

    file1=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/", tissue.names[t],
                "_AllPC/List.models.", tissue.names[t], ".all.RData", sep = "")

    file2=paste("/mnt/ceph/jarredk/AddisReRunFiles/List.models.",
                tissue.names[t], ".all.RData", sep = "")

    file3=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",
                tissue.names[t],"_AllPC/data.snp.cis.trans.final.",
                tissue.names[t],".V8.unique.snps.RData", sep = "")



    list.Lond=loadRData(fileName=file1)
    list.Addis=loadRData(fileName=file2)
    trios=loadRData(fileName=file3)


    ss=dim(trios)[2]/3

    trio.cn=as.data.frame(matrix(colnames(trios), nrow = ss, ncol = 3, byrow = T))

    tiss.table=as.data.frame(matrix(NA, nrow = ss, ncol = 16))
    colnames(tiss.table)=c("SNP", "cis.gene", "cis.gene.name", "trans.gene", "trans.gene.name",
                           "cis.chr", "trans.chr","cis.gene.start", "cis.gene.end", "trans.gene.start",
                           "trans.gene.end", "cis.gene.type", "trans.gene.type", "model", "Tissue", "Trio.Index")

    #fill in easy fields of master table
    tiss.table$SNP=trio.cn[,1]
    tiss.table$cis.gene=trio.cn[,2]
    tiss.table$trans.gene=trio.cn[,3]
    tiss.table$Tissue=rep(tissue.names[t], ss)
    tiss.table$Trio.Index=c(1:ss)

    if(which.method=="LOND"){
      models.list=list.Lond
    }else{
      models.list=list.Addis
    }


    tiss.table$model=find.models(models.list, vec.size = ss)

    cis.table=bm.data[match(tiss.table$cis.gene, bm.data$Stable.ID),]
    row.names(cis.table)=c(1:ss)
    trans.table=bm.data[match(tiss.table$trans.gene, bm.data$Stable.ID),]
    row.names(trans.table)=c(1:ss)

    tiss.table.filled=fill.table(tiss.table, cis.table, trans.table, size = ss)

    # print(cis.table[1:5,])
    # print(trans.table[1:5,])
    # print(head(tiss.table.filled))
    master.list[[t]]=tiss.table.filled
  }


  close(lbar2)
  print("...finished.")


  master.table=do.call(rbind.data.frame, master.list)

  if(save.table==TRUE){

    save(master.list, file=paste0("/mnt/ceph/jarredk/Manuscript/All_Trios_Master_table_",which.method,".RData"))


    write.table(master.table, file=paste0("/mnt/ceph/jarredk/Manuscript/All_Trios_Master_table_",which.method,".txt"),
                col.names = T, row.names = F, quote = F, sep = "\t")
  }

  return(master.table)

}

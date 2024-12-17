

#=================================================Helper_Function_1========================================================
get.unique.genesfile=function(input_tissue=NULL, save.file=TRUE){
  egenes.file.V8.1 <- read.delim (gzfile (paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/GTEx_Analysis_v8_eQTL/",
                                              input_tissue,".v8.egenes.txt.gz", sep = "")), header=TRUE)
  #dim(egenes.file.V8.1)


  egenes.file.V8 <- egenes.file.V8.1
  egenes.file.V8$gene_id <- gsub("\\..*","",egenes.file.V8.1$gene_id)
  egenes.file.V8$gene_name <- gsub("\\..*","",egenes.file.V8.1$gene_name)



  egenes.file.accending.V8 <- egenes.file.V8[order(egenes.file.V8['pval_nominal']), ]
  all.Genes.V8 <- egenes.file.accending.V8[,c(1,2,4,5,12,14,15,19,29,30)]
  dim(all.Genes.V8)
  all.Genes.V8[1:10,]



  list.of.eGenes.V8 <- which(all.Genes.V8[,9]<=0.05)
  eGenes.V8 <- all.Genes.V8[list.of.eGenes.V8,]
  dim(eGenes.V8)
  eGenes.V8[1:10,]

  #unique cis-snp
  eGenes.V8.unique <- eGenes.V8[!duplicated(eGenes.V8$variant_id), ]

  if(save.file==TRUE){
    save(eGenes.V8.unique,file=paste("/mnt/ceph/jarredk/AddisReRunFiles/",input_tissue,".eGenes.V8.unique.RData", sep = ""))
    
  }

  return(eGenes.V8.unique)

}








#================================================Helper_Function_2============================================

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}



#================================================START--HERE================================================
rerun.ADDIS=function(tissues.vec=NULL,tissues.vec2=NULL, save.files=TRUE, save.output=TRUE){

  
  library('MRPC', lib='/mnt/ceph/jarredk/Rpackages')

  #allocate for summarizing data structure 
  Total.Models.All.Tissues=as.data.frame(matrix(0, nrow = length(tissues.vec), ncol = 11))
  colnames(Total.Models.All.Tissues)=c("M0","M1","M2","M3","M4","Others","snps","genes","trios","samples","no.PCs")

  startingbar=rep("_", length(tissues.vec)*2)

  print(paste(startingbar, sep="",collapse=''))
  loading.bar=list()
  
  for(t in 1:length(tissues.vec)){
    #get .eGenes.V8.unique.Rdata file 
    eGenes.V8.unique=get.unique.genesfile(tissues.vec2[t], save.file = save.files)

    #load additional files
    file1=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",
          tissues.vec[t],"_AllPC/data.snp.cis.trans.final.",
          tissues.vec[t],".V8.unique.snps.RData", sep = "")

    file2=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",
                tissues.vec[t],"_AllPC/PCs.matrix.",
                tissues.vec[t], ".RData", sep="")

    file3=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissues.vec[t],
                "_AllPC/List.significant.asso1.", tissues.vec[t], ".RData", sep="")

    file4=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/PEER_Files_V8/Peerdata.", 
                tissues.vec[t], ".V8.RData", sep="")


    data.snp.cis.trans.final.V8.unique.snps=loadRData(fileName=file1)
    PCs.matrix=loadRData(fileName=file2)
    List.significant.asso1=loadRData(fileName=file3)
    Peerdata.V8=loadRData(fileName=file4)

    #check environment
    #print(ls())


    #----begin analysis-----#
    #truth for M0
    #V1-->T1
    Truth.M0 <- MRPCtruth$M0
    Adj.M0<- as(Truth.M0,"matrix")
    #V1-->T2
    Adj.M01 <- matrix(0,nrow=3,ncol = 3)
    rownames(Adj.M01) <- colnames(Adj.M01) <- colnames(Adj.M0)
    Adj.M01[1,3] <- 1

    #truth for M1
    #V1-->T1-->T2
    Truth.M1 <- MRPCtruth$M1
    Adj.M1<- as(Truth.M1,"matrix")
    #V1-->T2-->T1
    Adj.M11 <- matrix(0,nrow=3,ncol = 3)
    rownames(Adj.M11) <- colnames(Adj.M11) <- colnames(Adj.M1)
    Adj.M11[1,3] <- 1
    Adj.M11[3,2] <- 1

    #truth for M2
    #V1-->T1<--T2
    Truth.M2 <- MRPCtruth$M2
    Adj.M2<- as(Truth.M2,"matrix")
    #V1-->T2<--T1
    Adj.M21 <- matrix(0,nrow=3,ncol = 3)
    rownames(Adj.M21) <-colnames(Adj.M21) <-colnames(Adj.M2)
    Adj.M21[1,3] <- 1
    Adj.M21[2,3] <- 1
    #truth for M3
    #V1-->T1, V1-->T2
    Truth.M3 <- MRPCtruth$M3
    Adj.M3 <- as(Truth.M3,"matrix")

    #truth for M4
    #V1-->T1, V1-->T2, T1--T2
    Truth.M4 <- MRPCtruth$M4
    Adj.M4 <- as(Truth.M4,"matrix")


    #Initial results
    M00.FDR <- c()
    M11.FDR <- c()
    M22.FDR <- c()
    M33.FDR <- c()
    M44.FDR <- c()
    List.Associated.PCs <- list()

    # columns number for trios
    start.col <- seq(1,dim(data.snp.cis.trans.final.V8.unique.snps)[2],3)
    end.col <- seq(3,dim(data.snp.cis.trans.final.V8.unique.snps)[2],3)

    trios <- dim(data.snp.cis.trans.final.V8.unique.snps)[2]/3
    #print(trios)

    Match.PCs.col <- list()
    List.Match.significant.trios <- list()

    for (ii in 1:trios) {
      # data
      data <- data.snp.cis.trans.final.V8.unique.snps[,start.col[ii]:end.col[ii]]
      # search in each PCs
      for (m1 in 1:(dim(PCs.matrix)[1]-1)) {
        Match.PCs.col[[m1]] <- intersect(colnames(data),colnames(data.snp.cis.trans.final.V8.unique.snps)[List.significant.asso1[[m1]]])
      }
      List.Match.significant <- which(sapply(Match.PCs.col, function(e) length(e)!=0))
      List.Match.significant.trios[[ii]] <- List.Match.significant
      if(length(List.Match.significant)!=0)
      {
        data.withPC <- cbind(data,PCs.matrix[,List.Match.significant])
        colnames(data.withPC)[4:ncol(data.withPC)] <- paste("PC",List.Match.significant,sep="")
      }
      else
      {
        data.withPC <- data
      }
      
      
      
      n <- nrow (data.withPC)
      V <- colnames(data.withPC)     # Column names
      
      # Classical correlation
      suffStat <- list(C = cor(data.withPC,use = "complete.obs"),
                      n = n)
      
      MRPC.fit.FDR <- MRPC(data.withPC,
                          suffStat,
                          GV = 1,
                          FDR = 0.05,
                          indepTest = 'gaussCItest',
                          labels = V,
                          FDRcontrol = "ADDIS",
                          verbose = FALSE)
      
      #plot(MRPC.fit_FDR)
      Adj.infe1 <- as( MRPC.fit.FDR@graph,"matrix")
      Adj.infe <- Adj.infe1[1:3,1:3] #only consider snp, cis, trans
      colnames(Adj.infe) <- rownames(Adj.infe) <- colnames(Adj.M01)
      
      if(identical(Adj.M0,Adj.infe) || identical(Adj.M01,Adj.infe)){
        M00.FDR[ii]<-ii
      }
      
      if(identical(Adj.M1,Adj.infe) || identical(Adj.M11,Adj.infe)){
        M11.FDR[ii]<-ii
      }
      if(identical(Adj.M2,Adj.infe) || identical(Adj.M21,Adj.infe)){
        M22.FDR[ii]<-ii
      }
      if(identical(Adj.M3,Adj.infe)){
        M33.FDR[ii]<-ii
      }
      if(identical(Adj.M4,Adj.infe)){
        M44.FDR[ii]<-ii
      }
      
    }



  





  List.models.all <- list()
  List.models.all$M0 <- M00.FDR[!is.na(M00.FDR)]
  List.models.all$M1 <- M11.FDR[!is.na(M11.FDR)]
  List.models.all$M2 <- M22.FDR[!is.na(M22.FDR)]
  List.models.all$M3 <- M33.FDR[!is.na(M33.FDR)]
  List.models.all$M4 <- M44.FDR[!is.na(M44.FDR)]



  All.models.all<- matrix(c(length(List.models.all$M0),
                                length(List.models.all$M1),
                                length(List.models.all$M2),
                                length(List.models.all$M3),
                                length(List.models.all$M4),
                                length(setdiff(1:trios,c(List.models.all$M0,List.models.all$M1,
                                                          List.models.all$M2,List.models.all$M3,
                                                          List.models.all$M4))),
                                length(eGenes.V8.unique$variant_id),
                                dim(Peerdata.V8)[1],
                                trios,dim(Peerdata.V8)[2],
                                dim(PCs.matrix)[2]-1),
                              nrow = 1,ncol = 11)
  colnames(All.models.all) <- c("M0","M1","M2","M3","M4","Others","snps","genes","trios","samples","no.PCs")

  #print(All.models.all)


  Total.Models.All.Tissues[t,]=All.models.all
  #print(List.models.all)

  if(save.files==TRUE){
    
    save(All.models.all ,file=paste("/mnt/ceph/jarredk/AddisReRunFiles/All.models.",tissues.vec[t],".all.RData", sep = ""))
    save(List.models.all,file=paste("/mnt/ceph/jarredk/AddisReRunFiles/List.models.",tissues.vec[t],".all.RData",sep = ""))
    save(List.Match.significant.trios,file=paste("/mnt/ceph/jarredk/AddisReRunFiles/List.Match.significant.trios",t,".RData",sep = ""))
    
  }

  # #List.Match.significant.trios.1 <- which(sapply(List.Match.significant.trios, function(e) length(e)!=0))
  # PCasso.pc.all.trios <- which(sapply(List.Match.significant.trios, function(e) length(e)!=0))
  # length(PCasso.pc.all.trios)
  # #summary of associated PCs
  # Match.significant.trios.summary <- lengths(List.Match.significant.trios, use.names = TRUE)
  # summary(Match.significant.trios.summary)
  # #Total no of pc associated
  # List.significant.asso2 <- which(sapply(List.significant.asso1, function(e) length(e)!=0))
  # length(List.significant.asso2)



  loading.bar[[t]]=c(">",rep("=", t*2), startingbar[1:(length(startingbar)-t*2)])
  print(paste(paste(loading.bar[[t]], sep = "", collapse = ""), paste("...",round(((t*2)/(length(tissues.vec)*2))*100,2),"%", sep = ""), sep = ""))

  }

  if(save.output==TRUE){
  save(Total.Models.All.Tissues, file = "/mnt/ceph/jarredk/AddisReRunFiles/Total.Models.All.Tissues.Rdata")
  }

  print("done")


  return(Total.Models.All.Tissues)

    
  }
    
  #=========================================useful-lines=============================================== 
    #practice line of code
  #rerun.ADDIS(tissues.vec=c("WholeBlood","Vagina","Uterus"),tissues.vec2=c("Whole_Blood","Vagina","Uterus") ,save.files=FALSE) 
    
    
  #tissue.names=read.csv("/mnt/ceph/jarredk/AddisReRunFiles/tissuenames.csv", header = T)
  #output=rerun.ADDIS(tissues.vec=tissue.names$tissue.name1, tissues.vec2=tissue.names$tissue.name2, save.files=TRUE) 
    
    
    
  
    
    
  #output=rerun.ADDIS(tissues.vec=tissue.names[41,2], tissues.vec2=tissue.names[41,3], save.files=FALSE)
  


#=========================================================================================================
#-------------------------------------calculate-the-correlations------------------------------------------
#=========================================================================================================
#function to compute correlations between SNPS based on locations within each chromosome
#load in the necessary information:
#-------------------------------prep-data----------------------------------
loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#load in previously cleaned data
load(file="/mnt/ceph/jarredk/Methyl/Wksp1.Rdata")
genos.mat=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/fakegenosBIG.Rdata")
#EM.triolist=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/completedTrios1.Rdata")
#triobuildlist=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/triobuildlist.Rdata")
simulated.meta=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/fakegenoMeta.Rdata")
Mprobenames.final=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/Used.Mprobes.final.Rdata")
#load in methylation meta-data
#meta_M=read.csv(file = "/mnt/ceph/megheib/M_G_data/GPL13534_M.csv")
meta_M=read.csv(file = "/mnt/ceph/jarredk/Methyl/GPL13534_HumanMethylation450_15017482_v.1.1.csv")
#load in expression meta-data from biomart
meta_E=read.table(file = "/mnt/ceph/jarredk/Methyl/mart_export_GB37.txt", sep = ",", header = T)

#retrieve relevant metadata for the "In-Use" methylation probes:
imp.meta_M=na.omit(cbind(meta_M[,c(1,12)], as.numeric(meta_M[,13])))
matched.meta_M=match(Mprobenames.final, imp.meta_M[,1])
imp.meta_M.final=imp.meta_M[matched.meta_M,]
#retrieve relevant metadata for the "In-Use" Expression probes:
matched.meta_E=match(express.genenames, meta_E[,1])
imp.meta_E.final=meta_E[na.omit(matched.meta_E), ]
#filter out columns (genes) that didn't have chr position information according to BioMart
testexpress=express.aligned[,-attr(na.omit(matched.meta_E), "na.action")]
testexpress.gn=express.genenames[-attr(na.omit(matched.meta_E), "na.action")]
colnames(testexpress)=testexpress.gn
#-------------------------------------------------------------------------
#=========================helper_function_1===============================

cc.loop=function(probe.matched=NULL, Geno.matched=NULL, probe.mat=NULL, SNP.mat=NULL, bp.range=1000000, fn=NULL, verbose=FALSE ){
  
  for(j in 2:dim(probe.matched)[1]){
    
    #index the snps that are close to each Methylation probe
    snps.in.rangeM=which(abs(probe.matched[j,3]-Geno.matched[,2])<bp.range)
    Mprobe=rep(colnames(probe.mat)[j], length(snps.in.rangeM))
    SNP=colnames(SNP.mat)[snps.in.rangeM]
    cors=as.vector(cor(probe.mat[,j], SNP.mat[, snps.in.rangeM], use = "complete.obs"))
    add.tab=cbind.data.frame(Mprobe, SNP, cors)
    
    write.table(add.tab, file = fn,
                row.names = FALSE,
                col.names = FALSE,
                append = TRUE)
    
    if(isTRUE(verbose)){print(head(add.tab))}
  }
}

#=========================primary_function_1==============================

calc.corsV2=function(mmat=NULL, emat=NULL, gmat=NULL, GMInfo=NULL, GEInfo=NULL, genoInfo=NULL, chrs=NULL, bp.range=1000000){
  #SYNTAX:
  
  #GMInfo -- Gene Methylation matrix additional info in the form
  # [ probe_name  Chr  start_position  ]   
  
  #GEInfo -- Gene Expression matrix of additional info in the form
  # [ Gene_name  Chr  Start_position  ]
  
  #genomat -- the matrix/df of genotype data (with columnnames)
  
  #genoInfo -- SNP matrix/df of additional info in the form
  # [ Chr  Start_position  ]
  
  #Chrs -- passed from calc.cors
  
  for(i in 1:length(chrs)){
    
    snps.in.range=NULL
    snps.in.range2=NULL
    
    #get the idx's for each mat for chr1 
    Gidx=which(genoInfo[,1]==chrs[i])
    Midx=which(GMInfo[,2]==chrs[i])
    Eidx=which(GEInfo[,2]==chrs[i])
    
    #reduce meta data matrices
    chr.geno.matched=genoInfo[Gidx,]
    chr.GM.matched=GMInfo[Midx,]
    chr.GE.matched=GEInfo[Eidx, ]
    #print(chr.GE.matched)
    
    #reduce data matrices
    SNP.mat=gmat[,Gidx]
    methyl.mat=mmat[,Midx]
    express.mat=emat[,Eidx]
    
    #-----------------------Methylation-Probes-------------------------
    
    snps.in.rangeM=which(abs(chr.GM.matched[1,3]-chr.geno.matched[,2])<bp.range)
    Mprobe=rep(colnames(methyl.mat)[1], length(snps.in.rangeM))
    SNP=colnames(SNP.mat)[snps.in.rangeM]
    cors=as.vector(cor(methyl.mat[,1], SNP.mat[, snps.in.rangeM], use = "complete.obs"))
    init.table=cbind.data.frame(Mprobe, SNP, cors)
    colnames(init.table)=c("Methyl_Probe_ID", "SNP_ID", "Cor")
    
    print(dim(init.table))
    
    #fileName
    fnM=paste("/mnt/ceph/jarredk/Methyl/cor_Lists_and_Tables/", 
             paste0("M_chr_",i), "correlations.txt", sep = "")
    
    #initiate table
    write.table(init.table, file = fnM,
                row.names = FALSE)
    
    cc.loop(probe.matched=chr.GM.matched, 
            Geno.matched=chr.geno.matched, 
            probe.mat=methyl.mat, 
            SNP.mat=SNP.mat, 
            bp.range=bp.range, 
            fn=fnM, 
            verbose=FALSE)
    
    
    #--------------------------Expression-Probes----------------------------
    
    snps.in.rangeE=which(abs(chr.GE.matched[1,3]-chr.geno.matched[,2])<bp.range)
    Eprobe=rep(colnames(express.mat)[1], length(snps.in.rangeE))
    SNP=colnames(SNP.mat)[snps.in.rangeE]
    cors2=as.vector(cor(express.mat[,1], SNP.mat[, snps.in.rangeE], use = "complete.obs"))
    init.table2=cbind.data.frame(Eprobe, SNP, cors2)
    colnames(init.table2)=c("Express_Probe_ID", "SNP_ID", "Cor")
    
    print(dim(init.table2))
    
    #fileName
    fnE=paste("/mnt/ceph/jarredk/Methyl/cor_Lists_and_Tables/", 
              paste0("E_chr_",i), "correlations.txt", sep = "")
    
    #initiate table
    write.table(init.table2, 
                file = fnE,
                row.names = FALSE)
    
    cc.loop(probe.matched=chr.GE.matched, 
            Geno.matched=chr.geno.matched, 
            probe.mat=express.mat, 
            SNP.mat=SNP.mat, 
            bp.range=bp.range, 
            fn=fnE, 
            verbose=FALSE )
    
  }
  
}

start.time=Sys.time()
trycors1=calc.corsV2(mmat = methyl.resids2,
                   emat = testexpress,
                   gmat = genos.mat, 
                   GMInfo = imp.meta_M.final, 
                   GEInfo = imp.meta_E.final, 
                   genoInfo = simulated.meta, 
                   chrs=c("1"),
                   bp.range = 1000000)
end.time=Sys.time()

print(end.time-start.time)

















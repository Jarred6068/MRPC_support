


#=========================================================================================================
#-------------------------------------calculate-the-correlations------------------------------------------
#=========================================================================================================
#function to compute correlations between SNPS based on locations within each chromosome
#load in the necessary files needed to run function:
#-------------------------------prep-data----------------------------------
loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#load in previously cleaned data
genotype_data=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/fakegenosBIG.Rdata")
simulated.Geno_meta=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/fakegenoMeta.Rdata")
M_metadata=loadRData(imp.meta_M.final, file = "/mnt/ceph/jarredk/Methyl/Correlation_Code/meta_M.final.Rdata")
E_metadata=loadRData(imp.meta_E.final, file = "/mnt/ceph/jarredk/Methyl/Correlation_Code/meta_E.final.Rdata")
Edata=loadRData(testexpress, file = "/mnt/ceph/jarredk/Methyl/Correlation_Code/Expression_data_BM_aligned.final.Rdata")
Mdata=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/MethylData.RegressResids.Rdata")
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
  #mmat - the methylation matrix/df with column-names as probe ID's
  
  #emat - the expression matrix (aligned with BIOmart meta-data in GEInfo)
  
  #gmat -- the matrix/df of genotype data (with column-names as SNP ID's)
  
  #GMInfo -- Gene Methylation matrix additional info in the form
  # [ probe_name  Chr  start_position  ]   
  
  #GEInfo -- Gene Expression matrix of additional info in the form
  # [ Gene_name  Chr  Start_position  ]
  
  #genomat -- the matrix/df of genotype data (with columnnames)
  
  #genoInfo -- SNP matrix/df of additional info in the form
  # [ Chr  Start_position  ]
  
  #Chrs -- the chromosome(s) for which the script should be run can be
  #        can be a single or multi element vector
  
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
    
    #loop which writes correlations to above table
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
    
    #loop which writes correlations to above table
    cc.loop(probe.matched=chr.GE.matched, 
            Geno.matched=chr.geno.matched, 
            probe.mat=express.mat, 
            SNP.mat=SNP.mat, 
            bp.range=bp.range, 
            fn=fnE, 
            verbose=FALSE )
    
  }
  
}


#Example Run
start.time=Sys.time()
trycors1=calc.corsV2(mmat = Mdata,
                     emat = Edata,
                     gmat = genotype_data, 
                     GMInfo = M_metadata, 
                     GEInfo = E_metadata, 
                     genoInfo = simulated.meta, 
                     chrs=c("1"),
                     bp.range = 1000000)
end.time=Sys.time()

print(end.time-start.time)



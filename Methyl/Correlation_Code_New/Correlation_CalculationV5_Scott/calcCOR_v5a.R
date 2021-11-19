


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

#-------------------------------------------------------------------------
#=========================helper_function_1===============================

cc.loop=function(probe.matched=NULL, Geno.matched=NULL, probe.mat=NULL, SNP.mat=NULL, bp.range=1000000, fn=NULL){
  
  for(j in 2:dim(probe.matched)[1]){
    
    #index the snps that are close to each Methylation probe
    snps.in.rangeM=which(abs(probe.matched$coordinate[j]-Geno.matched$coordinate)<bp.range)
    
    if(length(snps.in.rangeM)==0){
      
      add.tab=c(colnames(probe.mat)[j], rep(NA,6))
      
      write.table(t(add.tab), file = fn,
                  row.names = FALSE,
                  col.names = FALSE,
                  append = TRUE,
                  sep = "\t",
                  quote = FALSE)
      
      
      
    }else{
      
      probe=rep(colnames(probe.mat)[j], length(snps.in.rangeM))
      SNP=colnames(SNP.mat)[snps.in.rangeM]
      cors=round(as.vector(cor(probe.mat[,j], SNP.mat[,snps.in.rangeM], use = "complete.obs")), 6)
      add.tab=cbind.data.frame(probe, 
                               SNP, 
                               cors,
                               Geno.matched$chr[snps.in.rangeM],
                               Geno.matched$coordinate[snps.in.rangeM],
                               Geno.matched$CountedAllele[snps.in.rangeM],
                               Geno.matched$OtherAllele[snps.in.rangeM])
      
      
      
      
      
      write.table(add.tab, file = fn,
                  row.names = FALSE,
                  col.names = FALSE,
                  append = TRUE,
                  sep = "\t",
                  quote = FALSE)
    }
    
  }
}

#=========================primary_function_1==============================

calc.cors=function(mmat=NULL, emat=NULL, gmat=NULL, GMInfo=NULL, GEInfo=NULL, genoInfo=NULL, chrs=NULL, fn=NULL, bp.range=1000000,
                   snps.as.rows=TRUE){
  #SYNTAX:
  #mmat - the methylation matrix/df with column-names as probe ID's
  
  #emat - the expression matrix (aligned with BIOmart meta-data in GEInfo)
  
  #gmat -- the matrix/df of genotype data: see README for formatting 
  
  #GMInfo -- Gene Methylation matrix additional info in the form
  # [ probe_name  Chr  start_position  ]   
  
  #GEInfo -- Gene Expression matrix of additional info in the form
  # [ Gene_name  Chr  Start_position  ]
  
  #genomat -- the matrix/df of genotype data (with columnnames)
  
  #genoInfo -- SNP matrix/df of additional info in the form
  # [ Chr  Start_position  ]
  
  #Chrs -- the chromosome(s) for which the script should be run can be
  #        can be a single or multi element vector
  
  #align all rows of genotype data to methylation data
  V=strsplit(row.names(mmat), "/")
  VV=matrix(unlist(V), nrow=dim(mmat)[1], ncol=2, byrow=TRUE)
  
  print("Extracting SNPs on this chromosome...")
  if(snps.as.rows==TRUE){
    
    message("User specified SNPs as rows: first column of data should contain SNP IDs")
    
    if(isFALSE(is.character(gmat[,1]))){
      stop(paste0("The first column of gmat data did not contain SNP IDs! See README.txt"))
    }
    
    idx=match(VV[,1], colnames(gmat))
    snp_ids=gmat[,1]
    #gmat=gmat[,-1]
    gmat=gmat[,na.omit(idx)]
    
  }else{
    
    message("User specified SNPs as columns: first column of data should contain BSGS IDs")
    
    if(isFALSE(is.character(gmat[,1]))){
      stop(paste0("The first column of gmat did not contain BSGS IDs! See README.txt"))
    }
    
    idx=match(VV[,1], gmat[,1])
    snp_ids=colnames(gmat)[-1]
    gmat=gmat[na.omit(idx),]
    gmat=gmat[,-1]
    
  }

  
  
  for(i in 1:length(chrs)){
    
    #get the idx's for each mat for chr1 
    Gidx=which(genoInfo$chr == chrs[i])
    Midx=which(GMInfo$chr == chrs[i])
    Eidx=which(GEInfo$chr == chrs[i])
    
    
    
    #reduce meta data matrices
    chr.geno.matched=genoInfo[Gidx,]
    chr.GM.matched=GMInfo[Midx,]
    chr.GE.matched=GEInfo[Eidx,]
    
    #reduce data matrices
    if(snps.as.rows==TRUE){
      SNP.mat=t(gmat[na.omit(match(chr.geno.matched[,1], snp_ids[Gidx])),])
    }else{
      SNP.mat=gmat[,na.omit(match(chr.geno.matched[,1], snp_ids[Gidx]))]
    }
    
    colnames(SNP.mat)=snp_ids[Gidx]
    methyl.mat=mmat[,na.omit(match(chr.GM.matched[,1], colnames(mmat)))]
    express.mat=emat[,na.omit(match(chr.GE.matched[,1], colnames(emat)))]
    
    
    #warnings to prevent mismatches
    #geno data
    if(length(which(genoInfo$chrs == chrs[i]))>0)
    
      if(isFALSE(isTRUE(all.equal(chr.geno.matched[,1], row.names(SNP.mat))))){
        warning(paste0("SNP data wasn't aligned!:  ", 
                      "SNP matrix contains ", dim(SNP.mat)[1])," SNPs and  metadata contains ", dim(chr.geno.matched)[1])
        stop(paste0("SNP probe names do not match metadata ID's for chr ==", chrs[i] ,"!"))
      }
    
    #methyl data
    
      if(isFALSE(isTRUE(all.equal(chr.GM.matched[,1], colnames(methyl.mat))))){
        warning(paste0("Methylation data isn't aligned  ", 
                      "M matrix contains ", dim(methyl.mat)[2])," probes and metadata contains ", dim(chr.GM.matched)[1])
        stop(paste0("Methylation probe names do not match metadata ID's for chr ==", chrs[i] ,"!"))
      }
    
    
    #expression data
    
      if(isFALSE(isTRUE(all.equal(chr.GE.matched[,1], colnames(express.mat))))){
        warning(paste0("Expression data isn't aligned!:  ", 
                     "E matrix contains ", dim(express.mat)[2]," genes and metadata contains ", dim(chr.GE.matched)[1]))
        stop(paste0("Expression probe names do not match metadata ID's for chr ==", chrs[i] ,"!"))
      }
    
    
    #-----------------------Methylation-Probes-------------------------
    
    print("Calculating Correlations With Methylation Probes...")
    snps.in.rangeM=which(abs(chr.GM.matched$coordinate[1]-chr.geno.matched$coordinate)<bp.range)
    
    if(length(snps.in.rangeM)==0){
      
      init.table=c(colnames(methyl.mat)[1], rep(NA, 6))
      
      names(init.table)=c("Methyl_Probe_ID", "SNP_ID", "Cor", "Chr", "Coordinate", "CountedAllele", "OtherAllele")
      
      
      #fileName
      fnM=paste(fn, paste0("M_chr",chrs[i]),"_", "correlations.txt", sep = "")
      
      #initiate table
      write.table(t(init.table), file = fnM,
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE)
      
    }else{
      
      Mprobe=rep(colnames(methyl.mat)[1], length(snps.in.rangeM))
      SNP=colnames(SNP.mat)[snps.in.rangeM]
      cors=round(as.vector(cor(methyl.mat[,1], SNP.mat[,snps.in.rangeM], use = "complete.obs")),6)
      init.table=cbind.data.frame(Mprobe, 
                                  SNP, 
                                  cors, 
                                  chr.geno.matched$chr[snps.in.rangeM],
                                  chr.geno.matched$coordinate[snps.in.rangeM],
                                  chr.geno.matched$CountedAllele[snps.in.rangeM],
                                  chr.geno.matched$OtherAllele[snps.in.rangeM])
      
      
      
      
      colnames(init.table)=c("Methyl_Probe_ID", "SNP_ID", "Cor", "Chr", "Coordinate", "CountedAllele", "OtherAllele")
      
      
      #fileName
      fnM=paste(fn, paste0("M_chr",chrs[i]),"_", "correlations.txt", sep = "")
      
      #initiate table
      write.table(init.table, file = fnM,
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE)
    }
    
    #loop which writes correlations to above table
    cc.loop(probe.matched=chr.GM.matched, 
            Geno.matched=chr.geno.matched, 
            probe.mat=methyl.mat, 
            SNP.mat=SNP.mat, 
            bp.range=bp.range, 
            fn=fnM)
    
    #--------------------------Expression-Probes----------------------------
    
    print("Calculating Correlations With Expression Probes...")
    snps.in.rangeE=which(abs(chr.GE.matched$coordinate[1]-chr.geno.matched$coordinate)<bp.range)

    if(length(snps.in.rangeE)==0){

      init.table=c(colnames(express.mat)[1], rep(NA, 6))

      names(init.table)=c("Express_Gene_ID", "SNP_ID", "Cor", "Chr", "Coordinate", "CountedAllele", "OtherAllele")


      #fileName
      fnE=paste(fn, paste0("E_chr",chrs[i]),"_", "correlations.txt", sep = "")

      #initiate table
      write.table(t(init.table), file = fnE,
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE)

    }else{

      Eprobe=rep(colnames(express.mat)[1], length(snps.in.rangeE))
      SNP=colnames(SNP.mat)[snps.in.rangeE]
      cors2=round(as.vector(cor(express.mat[,1], SNP.mat[,snps.in.rangeE], use = "complete.obs")),6)
      init.table2=cbind.data.frame(Eprobe,
                                   SNP,
                                   cors2,
                                   chr.geno.matched$chr[snps.in.rangeE],
                                   chr.geno.matched$coordinate[snps.in.rangeE],
                                   chr.geno.matched$CountedAllele[snps.in.rangeE],
                                   chr.geno.matched$OtherAllele[snps.in.rangeE])



      colnames(init.table2)=c("Express_Gene_ID", "SNP_ID", "Cor", "Chr", "Coordinate", "CountedAllele", "OtherAllele")


      #fileName
      fnE=paste(fn, paste0("E_chr",chrs[i]),"_", "correlations.txt", sep = "")

      #initiate table
      write.table(init.table2,
                  file = fnE,
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE)

    }

    #loop which writes correlations to above table
    cc.loop(probe.matched=chr.GE.matched,
            Geno.matched=chr.geno.matched,
            probe.mat=express.mat,
            SNP.mat=SNP.mat,
            bp.range=bp.range,
            fn=fnE)
    
  }
  
}







loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

filepath = "/mnt/ceph/jarredk/Methyl/CC_all/Correlation_CalculationV5_Scott/"
input.location = paste (filepath, "Input/", sep='')

Ecc=read.table(file="/mnt/ceph/jarredk/Methyl/Correlation_Calculation_Data/E_chr21_correlations.txt", header = T, sep="\t")
colnames(Ecc)=c("Express_Gene_ID", "SNP_ID", "Cor", "Chr", "Coordinate", "CountedAllele", "OtherAllele")

Mcc=read.table(file="/mnt/ceph/jarredk/Methyl/Correlation_Calculation_Data/M_chr21_correlations.txt", header = T, sep="\t")

gb37=read.table(file="/mnt/ceph/jarredk/Methyl/mart_export_GB37.txt", header = T, sep = ",")
colnames(gb37) = c("Gene.name", "chr", "start.bp", "end.bp", "IlmnID")

M_metadata=loadRData(fileName = paste(input.location, "All.Mprobes.Metadata.Rdata", sep = ""))
# Expression metadata: a meta data matrix of xxx probes and 3 columns [gene_ID     chr       coordinate]
E_metadata=loadRData(fileName = paste(input.location, "Express.BM.aligned.Ilmn.aligned.MetaData.Rdata", sep = ""))

hm450_meta=read.csv(file="/mnt/ceph/jarredk/Methyl/GPL13534_HumanMethylation450_15017482_v.1.1.csv")




check.E=function(number.of.its=1000, ECCfile=Ecc, biomart=gb37, verbose=TRUE){
  
  TF.vec=NULL
  
  for(i in 1:number.of.its){
    
    idx=sample(c(1:dim(ECCfile)[1]),1)
    
    if(verbose==TRUE){
      
      print(ECCfile[idx,])
      
      print(gb37[which(gb37$IlmnID==ECCfile$Express_Gene_ID[idx]), ])
      
    }
    
               
    TF.vec[i]=isTRUE(gb37$chr[which(gb37$IlmnID==ECCfile$Express_Gene_ID[idx])][1]==ECCfile$Chr[idx])  
    
    
    
  }
  
  print(sum(TF.vec))
  
}



check.M=function(number.of.its=1000, MCCfile=Mcc, hm450=hm450_meta, verbose=TRUE){
  
  TF.vec=NULL
  
  for(i in 1:number.of.its){
    
    idx=sample(c(1:dim(MCCfile)[1]),1)
    
    if(verbose==TRUE){
      
      print(MCCfile[idx,])
      
      print(hm450[which(hm450$IlmnID==MCCfile$Methyl_Probe_ID[idx]), 1:10])
      
    }
    
    
    TF.vec[i]=isTRUE(hm450$CHR[which(hm450$IlmnID==MCCfile$Methyl_Probe_ID[idx])][1]==MCCfile$Chr[idx])  
    
    
    
  }
  
  print(sum(TF.vec))
  
}



source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")

library('qvalue', lib="/mnt/ceph/jarredk/Rpackages")

top5=c(1,6, 33, 40, 48)
path='/mnt/ceph/jarredk/Reg_Net/'
tissues.vec=tissue.names[top5, 2:3]




#function to run post processing on GMAC output

run.postpoc=function(output.pvals=NULL, trio.ref=NULL){
  
  qvals=qvalue(output.pvals, fdr.level = 0.1)
  
  num.inferred=sum(qvals$significant)
  
  sig.trios=trio.ref[qvals$significant,]
  
  result=cbind.data.frame(qvals$qvalue[qvals$significant], output.pvals[qvals$significant])
  
  
  return(list(total.inferred=num.inferred, sig.trios=sig.trios, pq.values=result))
}

#load in output for the cis trios:

output.WB=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/WholeBlood/all_trios_output_cis50.Rdata')
output.AS=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/AdiposeSubcutaneous/all_trios_output_cis50.Rdata')
output.ArT=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/ArteryTibial/all_trios_output_cis50.Rdata')
output.MS=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/MuscleSkeletal/all_trios_output_cis50.Rdata')
output.SSE=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/SkinSunExposed/all_trios_output_cis50.Rdata')

#post process GMAC-cis-mediated
result.WB=run.postpoc(output.pvals = output.WB$output.table[,5], trio.ref = output.WB$output.table[,1:3])
result.AS=run.postpoc(output.pvals = output.AS$output.table[,5], trio.ref = output.AS$output.table[,1:3])
result.ArT=run.postpoc(output.pvals = output.ArT$output.table[,5], trio.ref = output.ArT$output.table[,1:3])
result.MS=run.postpoc(output.pvals = output.MS$output.table[,5], trio.ref = output.MS$output.table[,1:3])
result.SSE=run.postpoc(output.pvals = output.SSE$output.table[,5], trio.ref = output.SSE$output.table[,1:3])


#read in output for the trans trios:

#load in output for the cis trios:

output.t.WB=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/WholeBlood/all_trios_output_trans50.Rdata')
output.t.AS=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/AdiposeSubcutaneous/all_trios_output_trans50.Rdata')
output.t.ArT=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/ArteryTibial/all_trios_output_trans50.Rdata')
output.t.MS=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/MuscleSkeletal/all_trios_output_trans50.Rdata')
output.t.SSE=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/SkinSunExposed/all_trios_output_trans50.Rdata')

result.t.WB=run.postpoc(output.pvals = output.t.WB$output.table[,5], trio.ref = output.t.WB$output.table[,1:3])
result.t.AS=run.postpoc(output.pvals = output.t.AS$output.table[,5], trio.ref = output.t.AS$output.table[,1:3])
result.t.ArT=run.postpoc(output.pvals = output.t.ArT$output.table[,5], trio.ref = output.t.ArT$output.table[,1:3])
result.t.MS=run.postpoc(output.pvals = output.t.MS$output.table[,5], trio.ref = output.t.MS$output.table[,1:3])
result.t.SSE=run.postpoc(output.pvals = output.t.SSE$output.table[,5], trio.ref = output.t.SSE$output.table[,1:3])


result.table=cbind.data.frame(tissues.vec[,1], 
                              c(dim(output.AS$output.table)[1],
                                dim(output.ArT$output.table)[1],
                                dim(output.MS$output.table)[1],
                                dim(output.SSE$output.table)[1],
                                dim(output.WB$output.table)[1]),
                              c(result.AS$total.inferred,
                                result.ArT$total.inferred,
                                result.MS$total.inferred,
                                result.SSE$total.inferred,
                                result.WB$total.inferred),
                              c(result.t.AS$total.inferred,
                                result.t.ArT$total.inferred,
                                result.t.MS$total.inferred,
                                result.t.SSE$total.inferred,
                                result.t.WB$total.inferred),
                              c(123, 90, 104, 122, 131),
                              c(130, 113, 103, 123, 113))


colnames(result.table)=c("Tissue", 
                         "Total.Num.Samples", 
                         "GMAC.Total.M1T1",
                         "GMAC.Total.M1T2",
                         "LOND.total.M1", 
                         "ADDIS.total.M1")



write.csv(result.table, file="/mnt/ceph/jarredk/GMACanalysis/final_result_table.csv")



match.trios=function(tissues=tissue.vec[,1], which.mrpc="ADDIS"){
  
  for(i in 1:length(tissues)){
    
    
    #read in GMAC results
    output.cis=loadRData(fileName=paste('/mnt/ceph/jarredk/GMACanalysis/', tissues[i], '/all_trios_output_cis50.Rdata'))
    output.trans=loadRData(fileName=paste('/mnt/ceph/jarredk/GMACanalysis/', tissues[i], '/all_trios_output_trans50.Rdata'))
    
    #process 
    postpoc.cis=run.postpoc(output.pvals = output.cis$output.table[,5], 
                            trio.ref = output.cis$output.table[,1:3])
    
    postpoc.trans=run.postpoc(output.pvals = output.cis$output.table[,5], 
                              trio.ref = output.cis$output.table[,1:3])
    
    if(which.mrpc=="ADDIS"){
      
      
      
      
    }else{
      
      
      
    }
    
    
    
    
    
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
}








































































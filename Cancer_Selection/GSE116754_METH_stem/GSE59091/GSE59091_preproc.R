
source("/mnt/ceph/jarredk/Cancer_Selection/normalize_data.R")
#Stem cell methylation data consisting of 195 samples
#illumina infinium beadchip 450k 
setwd("/mnt/ceph/jarredk/Cancer_Selection/GSE116754_METH_stem/GSE59091/")
#read in raw data and meta_data
stem.meth2=read.table("GSE59091_unmethylated_methylated.txt", sep=" ", header=T)
stem.meth1.meta=read.csv("GSE59091_meta_info.csv")
dim(stem.meth2)


probe.id=stem.meth2$ID_REF
stem.meth.id.rm=stem.meth2[,-1]
#extract unmethylated
unmethylated=stem.meth.id.rm[,seq(1,(dim(stem.meth.id.rm)[2]-2), 3)]
unmethylated[1:5,1:8]
#extract methylated
methylated=stem.meth.id.rm[,seq(2,(dim(stem.meth.id.rm)[2]-1), 3)]
methylated[1:5,1:8]
#extract pvalues
detection.pvals=stem.meth.id.rm[,seq(3,(dim(stem.meth.id.rm)[2]), 3)]
detection.pvals[1:5,1:8]

prop.methylated=methylated/(methylated+unmethylated)

colnames(prop.methylated)=paste0("sample_",1:dim(prop.methylated)[2])
row.names(prop.methylated)=probe.id

prop.methylated[1:5,1:5]

write.table(prop.methylated, file = "GSE59091_prop_methyl.txt", quote = F, sep = "\t")
write.table(detection.pvals, file = "GSE59091_detection_pvals.txt", quote = F, sep = "\t")
#normalize props

norm.data=function(X, factor1, factor2){
  
  data=cbind.data.frame(y.log=X, Sex=as.factor(factor1), 
                        Reprog.method=as.factor(factor2))
  
  model=lm(y.log~Sex*Reprog.method, na.action = na.exclude, data=data)
  
  return(resid(model))
  
}



#get residueals

normalize.data=function(meth.data=NULL, meta.data=NULL, plot4=TRUE, file.id="GSE59091_"){
  
  pseudo=apply(meth.data, 2, pseudocount)
  write.table(pseudo, file = paste0(file.id,"pseudo.count.txt"), sep = "\t", quote = F)
  log.norm=apply(pseudo, 2, LT)
  write.table(log.norm, file = paste0(file.id,"log_normal.txt"), sep = "\t", quote=F)
  print(log.norm[1:5,1:8])
  #apply linear function to rows 
  resids=apply(t(log.norm), 2, 
               norm.data,
               factor1=meta.data$Sex, 
               factor2=meta.data$Reprogram_Method)
  
  print(t(resids[1:5,1:5]))
  
  # colnames(resids)=colnames(meth.data)
  # row.names(resids)=row.names(meth.data)
  
  
  
  write.table(t(resids), file = paste0(file.id,"residuals.txt"), sep = "\t", quote=F)
  
  x.sample=sample(c(1:dim(meth.data)[2]), 4)
  y.sample=sample(c(1:dim(meth.data)[2]), 4)
  
  if(plot4==TRUE){
    
    pdf(paste0(file.id,"plot.pdf"))
    par(mfrow=c(2,2))
    
    for(i in 1:4){
      
      x=c(t(meth.data)[,x.sample[i]], resids[,x.sample[i]])
      y=c(t(meth.data)[,y.sample[i]], resids[,y.sample[i]])
      colvar=c(rep("non.norm", dim(meth.data)[2]), rep("norm", dim(meth.data)[2]))
      
      plot(x,y, type="p", col=as.factor(colvar), pch=21, bg=as.factor(colvar),
           xlab = colnames(meth.data)[x.sample[i]],
           ylab = colnames(meth.data)[y.sample[i]],
           main = "Before and After transformation")
      
    }
    
    par(mfrow=c(1,1))
    dev.off()
    
  }
}
  

normalize.data(meth.data = prop.methylated, meta.data = stem.meth1.meta)


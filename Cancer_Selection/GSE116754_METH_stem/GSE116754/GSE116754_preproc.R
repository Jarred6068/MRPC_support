
#source in helper functions

source("/mnt/ceph/jarredk/Cancer_Selection/normalize_data.R")
#reading in second stem cell data set

#set working directory
setwd("/mnt/ceph/jarredk/Cancer_Selection/GSE116754_METH_stem/GSE116754/")
#load in data and meta data files
meth.data=read.table("GSE116754_methylated_and_unmethylated_signal.txt", sep="\t", header =T)
meta.data=read.csv("GSE116754_meta_info.csv")

#remove white space and replace with NA
meta.data.na=apply(meta.data, 2, function(x) replace(x, x=="", NA))
meta.data.na=as.data.frame(meta.data.na)
write.table(meta.data.na, file = "GSE116754_STEM_meta_data_no_na.txt", quote = F, sep = "\t")

#processing
probe.id=meth.data$ID_REF
meth.id.rm=meth.data[,-1]
#extract unmethylated
unmethylated=meth.id.rm[,seq(1,(dim(meth.id.rm)[2]-2), 3)]
unmethylated[1:5,1:8]
#extract methylated
methylated=meth.id.rm[,seq(2,(dim(meth.id.rm)[2]-1), 3)]
methylated[1:5,1:8]
#extract pvalues
detection.pvals=meth.id.rm[,seq(3,(dim(meth.id.rm)[2]), 3)]
detection.pvals[1:5,1:8]

prop.methylated=methylated/(methylated+unmethylated)

colnames(prop.methylated)=paste0("SAMPLE_",c(1:dim(methylated)[2]))
row.names(prop.methylated)=probe.id

prop.methylated[1:5,1:5]

write.table(prop.methylated, file = "GSE116754_STEM_prop_methyl.txt", quote = F, sep = "\t")
write.table(detection.pvals, file = "GSE116754_STEM_detection_pvals.txt", quote = F, sep = "\t")
#normalize props



norm.data=function(X, factor1){
  
  data=cbind.data.frame(y.log=X, Sample.Type=as.factor(factor1))
  
  model=lm(y.log~Sample.Type, na.action = na.exclude, data=data)
  
  return(resid(model))
  
}

#get residueals

normalize.data=function(meth.data=NULL, meta.data=NULL,which.covars=NULL, plot4=TRUE, file.id="GSE116754_STEM_"){
  
  #handle pseudocounts
  pseudo=apply(meth.data, 2, pseudocount)
  write.table(pseudo, file = paste0(file.id,"pseudo.count.txt"), sep = "\t", quote = F)
  #log transform
  log.norm=apply(pseudo, 2, LT)
  write.table(log.norm, file = paste0(file.id,"log_normal.txt"), sep = "\t", quote=F)
  print(log.norm[1:5,1:8])
  #apply linear function to rows and get residuals
  resids=apply(t(log.norm), 2, 
               norm.data,
               factor1 = meta.data.na$Sample_Type)
  
  print(t(resids[1:5,1:5]))
  
  #colnames(resids)=colnames(meth.data)
  #row.names(resids)=row.names(meth.data)
  
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


normalize.data(meth.data = prop.methylated, meta.data = meta.data.na)






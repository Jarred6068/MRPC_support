#####################################################################

#reading normal human breast tissue files

source("/mnt/ceph/jarredk/Cancer_Selection/normalize_data.R")
#total of 121 healthy human breast tissue samples:
#read in signal intensities matrix (methylated sites count, unmethylated sites ct, and pvalues)
setwd("/mnt/ceph/jarredk/Cancer_Selection/GSE101961_normal_breast/")
mat.sig=read.table("GSE101961_Matrix_signal.txt", header=T, sep="\t")
meta.data=read.csv("GSE101961_meta_info.csv")

#parse colnames
names.list=strsplit(colnames(mat.sig), split=".", fixed = TRUE)
#get the pval columns:
methyl.sig.pcols=lapply(names.list, FUN = function(x)  any(x=="Pval"))
mat.sig.pval.only=mat.sig[,c(1, which(methyl.sig.pcols==TRUE))]
dim(mat.sig.pval.only) #dimension should be probes x samples + ID column

#which columns are methylated sig cts -- allocate to separate matrix
methyl.sig.tf=lapply(names.list, FUN = function(x)  any(x=="Methylated"))
mat.sig.methyl=mat.sig[,c(1, which(methyl.sig.tf==TRUE))]
dim(mat.sig.methyl) #dimension should be probes x samples + ID column

#which columns are methylated sig cts -- same as above but for unmethylated
unmethyl.sig.tf=lapply(names.list, FUN = function(x)  any(x=="Unmethylated"))
mat.sig.unmethyl=mat.sig[,c(1,which(unmethyl.sig.tf==TRUE))]
dim(mat.sig.unmethyl)

#calculate % of methylation at each probe site for each sample

final.methyl=cbind.data.frame(mat.sig.methyl[,-1]/(mat.sig.methyl[,-1]+mat.sig.unmethyl[,-1]))
colnames(final.methyl)=colnames(mat.sig.methyl)[-1]
row.names(final.methyl)=mat.sig[,1]
#save
write.table(final.methyl, file="GSE101961.unnorm.methyl.final.txt", quote=F, sep="\t")
write.table(mat.sig.pval.only, file="GSE101961.signal.pvals.txt", quote=F, sep="\t")


##normalize data

#normalize props

norm.data=function(X, factor1, factor2, factor3){

  data=cbind.data.frame(y.log=X, Age=factor1, Race=as.factor(factor2), BMI=factor3)

  model=lm(y.log~poly(Age,2)*poly(BMI,2)*Race, na.action = na.exclude, data=data)

  return(resid(model))

}



#get residueals

normalize.data=function(meth.data=NULL, meta.data=NULL, plot4=TRUE, file.id="GSE101961_normal_breast_"){

  #handle pseudocounts
  pseudo=apply(meth.data, 2, pseudocount)
  write.table(pseudo, file = paste0(file.id,"pseudo.count.txt"), sep = "\t", quote = F)
  #log transform
  log.norm=apply(pseudo, 2, LT)
  write.table(log.norm, file = paste0(file.id,"log_normal.txt"), sep = "\t", quote=F)
  print(log.norm[1:5,1:8])
  #apply linear function to rows and get residuals
  resids=apply(t(log.norm), 2, norm.data,
               factor1=meta.data$Age,
               factor2=meta.data$Race,
               factor3=meta.data$BMI)

  print(t(resids[1:5,1:5]))

  # colnames(resids)=colnames(meth.data)
  # row.names(resids)=row.names(meth.data)

  write.table(t(resids), file = paste0(file.id,"residuals.txt"), sep = "\t", quote=T)

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

normalize.data(meth.data = final.methyl, meta.data = meta.data)


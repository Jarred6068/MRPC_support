

rraddis<- read.csv("C:/Users/Bruin/Desktop/Research Assistantship/ADDIS_ReRun/ReRunADDIS_Summary.csv")
rrlond<- read.csv("C:/Users/Bruin/Desktop/Research Assistantship/ADDIS_ReRun/ReRunLOND_Summary.csv")
nam=c("M0","M1","M2","M3","M4","Other")
modcts=rraddis[,c(6:11)]
colnames(modcts)=nam
modcts2=rrlond[,c(6:11)]
colnames(modcts2)=nam
modpcts=rraddis[,c(12:17)]
colnames(modpcts)=nam
modpcts2=rrlond[,c(12:17)]
colnames(modpcts2)=nam

L=rep("LOND", dim(modcts2)[1])
A=rep("ADDIS", dim(modcts)[1])

AplusL=c(A,L)

#pairs(modcts)
#pairs(modcts2)

plot.list=list()

par(mfrow=c(2,3))
for(i in 1:6){
  hist(modcts[,i],main = paste("histogram of:", nam[i], sep = " "))
}
par(mfrow=c(1,1))




library(ggplot2)
for(i in 1:6){

  df=cbind.data.frame( FDR=as.factor(AplusL), M=c(modcts[,i], modcts2[,i]))
  
  plot.list[[i]]=ggplot(df, aes(x=M, color=FDR, fill=FDR)) +geom_density(alpha=0.3)+
    ggtitle(label = paste("Density plot for:", nam[i], sep = " "))
    
}


library(gridExtra)
grid.arrange(grobs=plot.list,nrows=3, ncol=2)

#change table
average.change.ct=colMeans(modcts-modcts2)
average.change.ct.sd=apply(modcts-modcts2, 2, sd)
average.change.pct=colMeans(modpcts-modpcts2)
average.change.pct.sd=apply(modpcts-modpcts2,2,sd)

change.table=rbind.data.frame(average.change.ct, average.change.ct.sd, 
                              average.change.pct, average.change.pct.sd)
colnames(change.table)=nam
row.names(change.table)=c("mean.change.ct", "SD.change.ct", "mean.change%", "SD.change%")




#Bootstrap CI's

Addis.cIs=as.data.frame(matrix(0, nrow = 6, ncol = 3))
colnames(Addis.cIs)=c("lower.limit", "mean","upper.limit")
row.names(Addis.cIs)=nam
B=1000
alpha=0.05


#calculate confidence interval for each model type avg count(ADDIS)
for(i in 1:6){
  bsdata=rep(0,B)
  
  for(j in 1:B){
    
    bsdata[j]=mean(sample(modcts[,i], dim(modcts)[1], replace = TRUE))
    
  }
  
  #hist(bsdata)
  SE=sd(bsdata)
  M=mean(bsdata)
  lower=M-qnorm(alpha/2)*(SE/sqrt(length(bsdata)))
  upper=M+qnorm(alpha/2)*(SE/sqrt(length(bsdata)))
  
  Addis.cIs[i,]=round(c(lower,M,upper),3)
}



Lond.cIs=as.data.frame(matrix(0, nrow = 6, ncol = 3))
colnames(Lond.cIs)=c("lower.limit", "mean","upper.limit")
row.names(Lond.cIs)=nam
#calculate confidence interval for each model type avg count(LOND)
for(i in 1:6){
  bsdata=rep(0,B)
  
  for(j in 1:B){
    
    bsdata[j]=mean(sample(modcts2[,i], dim(modcts2)[1], replace = TRUE))
    
  }
  
  #print(bsdata)
  #hist(bsdata, breaks = 15)
  SE=sd(bsdata)
  M=mean(bsdata)
  lower=M-qnorm(alpha/2)*(SE/sqrt(length(bsdata)))
  upper=M+qnorm(alpha/2)*(SE/sqrt(length(bsdata)))
  
  Lond.cIs[i,]=round(c(lower,M,upper),3)
}






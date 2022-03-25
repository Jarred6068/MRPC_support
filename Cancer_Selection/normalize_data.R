

#############################################################
#Normalize methylation data in same method as BSGS data

LT=function(X){
  
  Z=ifelse(is.na(X), NA, log(X/(1-X)))
  
  return(Z)
}

pseudocount=function(X){
  
  if(length(which(X==0))>0){ X[which(X==0)]=0.001 }
  if(length(which(X==1))>0){ X[which(X==1)]=0.999 }
  
  return(X)
  
}


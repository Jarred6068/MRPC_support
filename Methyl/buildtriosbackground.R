


load(file="/mnt/ceph/jarredk/Methyl/Wksp1.Rdata")

start.time=Sys.time()

build.trios=function(emat=NULL, mmat=NULL, tbl=NULL, enames=NULL){
  
  
  trio.full.list=vector("list", length = length(enames))
  names(trio.full.list)=enames
  checker=integer(0)
  
  for(i in 1:length(enames)){
    print(i)
    
    if(identical(checker, tbl[[i]])){
      
      trio.full.list[[i]]=NA
      
    }else{
      
      #print(mmat[,tbl[[i]] ])
      submatM=as.matrix(mmat[, tbl[[i]] ])
      trio.gl=vector("list", length = length(tbl[[i]]) )
      
      for(j in 1:length(tbl[[i]])){
        combmat=cbind.data.frame(emat[,i], submatM[, j])
        colnames(combmat)=c(paste0(enames[i],";",colnames(emat)[i]), colnames(submatM)[j])
        row.names(combmat)=row.names(submatM)
        #print(combmat)
        trio.gl[[j]]=combmat
      }
      
      trio.full.list[[i]]=trio.gl
      
    }
    
    
  }
  
  return(trio.full.list)
  
}

trios=build.trios(emat = express.aligned, mmat = methyl.resids2, tbl = triobuildlist, enames = express.genenames)

end.time=Sys.time()

print(end.time-start.time)

save(trios, file = "/mnt/ceph/jarredk/Methyl/completedTrios1.Rdata")


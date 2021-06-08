

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

CNA=loadRData(fileName="C:/Users/Bruin/Documents/GitHub/MRPC_support/Misc/data.CNA.peer.RData")
EM=loadRData(fileName="C:/Users/Bruin/Documents/GitHub/MRPC_support/Misc/data.exp.peer.fake.RData")
MM=loadRData(fileName="C:/Users/Bruin/Documents/GitHub/MRPC_support/Misc/sub.n.n.TCGA.meth.t.2.RData")

CNA.names=colnames(CNA)
express.names=colnames(EM)
methyl.names.list=strsplit(colnames(MM), split = ";")



qc=function(X, V){ X==V }

match.gn.parsed=function(egn, mgnp){
  
  indexlist1=vector("list", length(egn))
  vec1=c(1:length(mgnp))
  vec11=NULL
  
  for(i in 1:length(egn)){
    
    loglist1=lapply(mgnp, qc, egn[i])
    
    vec11=vec1[unlist(lapply(loglist1, any))]
    
    indexlist1[[i]]=vec11
    
  }
  
  return(indexlist1)
  
}


try1=match.gn.parsed(express.names, methyl.names.list)

#clean up matrices and align rows.
M.probenames=MM[1,]

rm1=match("Probe.ID1", row.names(MM))
rm2=match("Genomic_Coordinate1", row.names(MM))
MM=MM[-c(1,2,3,rm1,rm2),]
colnames(MM)=M.probenames

align.rows=match(row.names(MM),row.names(EM))
MM2=MM[-which(is.na(align.rows)),]
MM2=apply(MM2, 2, as.numeric)
MM2[1:5,1:5]
EM2=EM[na.omit(align.rows),]
EM2[1:5,1:5]

build.trios1=function(emat=NULL, mmat=NULL, tbl=NULL, enames=NULL){
  #emat -- the expression data matrix
  #mmat -- the methylation data matrix
  #tbl -- the triobuildlist (output from match.gn.parsed)
  #enames -- the gene names that correspond to the columns of the 
  #           expression matrix
  
  trio.full.list=vector("list", length = length(enames))
  names(trio.full.list)=enames
  checker=integer(0)
  
  for(i in 1:length(enames)){

    
    if(identical(checker, tbl[[i]])){
      
      trio.full.list[[i]]=NA
      
    }else{
      
      #print(mmat[,tbl[[i]] ])
      submatM=as.matrix(mmat[, tbl[[i]] ])
      colnames(submatM)=colnames(mmat)[tbl[[i]] ]
      print(head(submatM))
      trio.gl=vector("list", length = length(tbl[[i]]) )
      
      for(j in 1:length(tbl[[i]])){
        combmat=cbind.data.frame(emat[,i], submatM[, j])
        colnames(combmat)=c(paste0(enames[i],";",colnames(emat)[i]), colnames(submatM)[j])
        row.names(combmat)=row.names
        #print(combmat)
        trio.gl[[j]]=combmat
      }
      
      trio.full.list[[i]]=trio.gl
      
    }
    
    
  }
  
  return(trio.full.list)
  
}


bt=build.trios(emat=EM2, mmat=MM2, tbl=try1, enames=express.names)


build.trios2=function(gmat=NULL, emat=NULL, mmat=NULL, tbl=NULL, enames=NULL){
  
  
  bt1=build.trios1(emat = emat, mmat = mmat, tbl = tbl, enames = enames)
  
  
  
  
}

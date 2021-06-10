

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

 CNA=loadRData(fileName="C:/Users/Bruin/Documents/GitHub/MRPC_support/Misc/data.CNA.peer.RData")
 EM=loadRData(fileName="C:/Users/Bruin/Documents/GitHub/MRPC_support/Misc/data.exp.peer.fake.RData")
 MM=loadRData(fileName="C:/Users/Bruin/Documents/GitHub/MRPC_support/Misc/sub.n.n.TCGA.meth.t.2.RData")

# CNA=loadRData(fileName="/mnt/ceph/jarredk/Misc/data.CNA.peer.RData")
# EM=loadRData(fileName="/mnt/ceph/jarredk/Misc/data.exp.peer.fake.RData")
# MM=loadRData(fileName="/mnt/ceph/jarredk/Misc/sub.n.n.TCGA.meth.t.2.RData")

CNA.names=colnames(CNA)
express.names=colnames(EM)

#remove Methyl columns with all or nearly all NA's

idxc=NULL
for(i in 1:dim(MM)[2]){if(sum(is.na(MM[,i]))>25){idxc[i]=i}else{idxc[i]=NA}}

MM=MM[,-na.omit(idxc)]

methyl.names.list=strsplit(colnames(MM), split = ";")

#get the trio construction list

#====================================function=======================================

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

#===================================================================================


try1=match.gn.parsed(express.names, methyl.names.list)

#clean up matrices and align rows.
M.probenames=MM[1,]

rm1=match("Probe.ID1", row.names(MM))
rm2=match("Genomic_Coordinate1", row.names(MM))
MM=MM[-c(1,2,3,rm1,rm2),]
colnames(MM)=M.probenames

#make sure all matrices have same number of rows
align.rows=match(row.names(MM),row.names(EM))
MM2=MM[-which(is.na(align.rows)),]
MM2=apply(MM2, 2, as.numeric)
MM2[1:5,1:5]
EM2=EM[na.omit(align.rows),]
CNA2=CNA[na.omit(align.rows),]
EM2[1:5,1:5]
CNA2[1:5,1:5]


#first layer of trio constructions
#====================================function=======================================

build.trios=function(emat=NULL, mmat=NULL, tbl=NULL, enames=NULL){
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
      #print(head(submatM))
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



#====================================function=======================================

#adds in genotype data to construct a list of full trios
#run and view the output object to see the format

bt=build.trios(emat=EM2, mmat=MM2, tbl=try1, enames=express.names)


build.trios2=function(gmat=NULL, emat=NULL, mmat=NULL, tbl=NULL, enames=NULL){
  
  trio.list.final=vector("list", length = dim(gmat)[2])
  
  bt1=build.trios(emat = emat, mmat = mmat, tbl = tbl, enames = enames)
  
  #rm NaN's
  idx=NULL
  for(i in 1:length(bt1)){ if(any(is.na(bt1[[i]]))){ idx[i]=i}else{idx[i]=NA} }
  
  bt1=bt1[-na.omit(idx)]
  
  for(i in 1:dim(gmat)[2]){

    pairs1=list()

    for(j in 1:length(bt1)){

      pairs1[[j]]=lapply(bt1[[j]], cbind.data.frame, gmat[,i])

    }
    
    names(pairs1)=names(bt1)
  
    trio.list.final[[i]]=pairs1

  }
  
  names(trio.list.final)=colnames(gmat)
  
  return(trio.list.final)
  
}



start.time=Sys.time()
bt2=build.trios2(gmat = CNA2, emat=EM2, mmat=MM2, tbl=try1, enames=express.names)
end.time=Sys.time()

print(end.time-start.time)


















































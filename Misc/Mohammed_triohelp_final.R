

#set directory

fn="/mnt/ceph/jarredk/Misc/"

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


#load in real data
CNA=loadRData(fileName=paste(fn,"CNA.data.Rdata",sep = ""))
EM=loadRData(fileName=paste(fn,"EM.data.Rdata",sep = ""))
MM=loadRData(fileName=paste(fn,"MM.data.Rdata",sep=""))

#match up rows
matched.rows=match(row.names(MM), row.names(CNA))

#reduce rows to be consistent across all 3 data
CNA=CNA[na.omit(matched.rows),]
EM=EM[na.omit(matched.rows),]
MM=MM[-attr(na.omit(matched.rows), "na.action"),]

# #filter columns with large number of NA's
# idxc=NULL
# for(i in 1:dim(MM)[2]){if(sum(is.na(MM[,i]))>25){idxc[i]=i}else{idxc[i]=NA}}
# 
# #apply filter
# MM=MM[,-na.omit(idxc)]

save(MM, file = paste(fn,"MM.data.formatted.Rdata",sep = ""))
save(EM, file = paste(fn,"EM.data.formatted.Rdata",sep = ""))
save(CNA, file = paste(fn,"MM.data.formatted.Rdata",sep = ""))

#get a list of the methylation gene/probe names
methyl.names.list=strsplit(colnames(MM), split = ";")
express.names=colnames(EM)

#get the trio construction list

#====================================function=======================================

qc=function(X, V){ X==V }

match.gn.parsed=function(egn, mgnp){
  
  indexlist1=vector("list", length(egn))
  vec1=c(1:length(mgnp))
  vec11=NULL
  
  for(i in 1:length(egn)){
    
    loglist1=lapply(mgnp, qc, egn[i])
    
    vec11=vec1[na.omit(unlist(lapply(loglist1, any)))]
    
    indexlist1[[i]]=vec11
    #print(indexlist1[[i]])
    #print(i)
  }
  
  return(indexlist1)
  
}

#===================================================================================

#----------------------formatting example data--------------------------------
try1=match.gn.parsed(express.names, methyl.names.list)

save(try1, file=paste(fn,"TBL.Rdata",sep = ""))

TBL=loadRData(fileName = paste(fn,"TBL.Rdata",sep = ""))

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
        colnames(combmat)=c(enames[i], colnames(submatM)[j])
        row.names(combmat)=row.names(submatM)
        #print(combmat)
        trio.gl[[j]]=combmat
        names(trio.gl)=colnames(submatM)[j]
      }
      
      names(trio.gl)=
        trio.full.list[[i]]=trio.gl
      
    }
    
    
  }
  
  return(trio.full.list)
  
}


bt=build.trios(emat=EM, mmat=MM, tbl=TBL, enames=express.names)

save(bt, file = paste(fn,"builttrios1.Rdata",sep = ""))
bt=loadRData(fileName = paste(fn, "builttrios1.Rdata",sep=""))



#this function constructs all possible trios for a given CNA probe and returns the shortened list of trios
#for this probe (removes NA elements in the list)

trios.comb=function(gmat=NULL, emat=NULL, mmat=NULL, bt1=NULL, enames=NULL){
  
  #   #rm NaN's from the list
  # idx=NULL
  # for(i in 1:length(bt1)){ if(any(is.na(bt1[[i]]))){ idx[i]=NA}else{idx[i]=NA} }
  # bt1=bt1[-na.omit(idx)]
  
  #pre allocate space
  list1=list()
  list2=list()
  
  for(i in 1:length(bt1)){

    for(j in 1:length(bt1[[i]])){
    
      #add in the probe to each pair of M and E data
    
      comb=cbind.data.frame(gmat[,i], bt1[[i]][[j]])
      colnames(comb)=c(colnames(gmat)[i], colnames(bt1[[i]][[j]]))
      #print(comb)
    
      list1[[j]]=comb
    
     }
   
    names(list1)=names(bt1[[i]])
    list2[[i]]=list1
  
   }
  
  #get the names of each trio
  namelist=list()
  for(i in 1:length(list2)){namelist[[i]]=lapply(list2[[i]], colnames)}
  
  #return both the list and it formated as a matrix
  names(list2)=names(bt1)
  
  V1=lapply(list2, as.data.frame)
  V2=as.data.frame(V1)
  #apply the namelist to cols
  colnames(V2)=unlist(namelist)
  
  #every 3 cols starting at mat[,1] is a trio
  return(list(list=V1, mat=V2))
  #return(list2)
  
}



#Ex
#check the format of the inputs to see how the matrices need to be arranged

test1=trios.comb(gmat=CNA, emat = EM, mmat = MM, bt1 = bt, enames = express.names)







































































































































































































































































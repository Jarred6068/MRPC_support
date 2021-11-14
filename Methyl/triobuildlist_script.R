


#A simple helper function for easier reading in of .Rdata files with a specific name
#================================================Helper_Function============================================

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#---------------------------------------------------------------------------------------------

E.genenames=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/ExpressData/Peerdata.bsgs_GeneIDs.Rdata")
UCSC_ref_gn_parsed=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/UCSC_ref_gn_parsed.Rdata")

#a function to match gene names from expression to the parsed list of methylation probe genes
#to pre-index the trio construction:
#this function is slow and can take some time (>10 hrs) to run

qc=function(X, V){ X==V }

match.gn.parsed=function(egn, mgnp){
  
  indexlist1=vector("list", length(egn))
  vec1=c(1:length(mgnp))
  vec11=NULL
  
  for(i in 1:length(egn)){
    
    loglist1=lapply(mgnp, qc, egn[i])
    
    vec11=vec1[unlist(lapply(loglist1, any))]
    
    indexlist1[[i]]=vec11
    print(indexlist1[[i]])
    
  }
  
  return(indexlist1)
  
}

#a list of the methylation probes matched to each gene in express.genenames
triobuildlist=match.gn.parsed(E.genenames, UCSC_ref_gn_parsed)

save(triobuildlist, file = "/mnt/ceph/jarredk/Methyl/ExpressData/triobuildlist.Rdata")


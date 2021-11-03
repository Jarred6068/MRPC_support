


#A simple helper function for easier reading in of .Rdata files with a specific name
#================================================Helper_Function============================================

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#---------------------------------------------------------------------------------------------
#READING IN DATA AND FORMATTING:

#read in normalized, filtered Expression Data
E=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/ExpressData/Peerdata.bsgs.RData")
#read in the associated Gene Names for the expression data
E.genenames=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/ExpressData/Peerdata.bsgs_GeneIDs.Rdata")
#read in normalized, filtered Methylation Data
M=loadRData(fileName="/mnt/ceph/jarredk/Methyl/MethylData.RegressResids.Rdata")
#read in the UCSC ref genenames for Methylation Data
E.genenames=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/MethylData/MethylData.RegressResids.UCSC_RefGene_Name.Rdata")
#read in BioMart Data containing Illumina ID's
gb37=read.table(file="/mnt/ceph/jarredk/Methyl/mart_export_GB37.txt", header = T, sep = ",")
#rename cols
colnames(gb37) = c("Gene.name", "chr", "start.bp", "end.bp", "IlmnID")
#Read in Human Methylation 450 data
hm450_meta=read.csv(file="/mnt/ceph/jarredk/Methyl/GPL13534_HumanMethylation450_15017482_v.1.1.csv")
#since the col dim of M < row dim of hm450, reduce hm450 according the probes actually present in M

hm450=hm450_meta[match(colnames(M), hm450_meta$IlmnID),]

#check to make sure dims match
dim(hm450)
dim(M)

#replace "" in hm450$UCSC_RefGene_Name with NA
ucsc_ref=hm450$UCSC_RefGene_Name
for(i in 1:dim(hm450)[1]){
  
  ucsc_ref[i]=ifelse(ucsc_ref[i]=="", NA, hm450$UCSC_RefGene_Name[i])
  
}

print(ucsc_ref[1:10])
print(hm450$UCSC_RefGene_Name[1:10])

#reduce hm450 by removing rows with missing UCSC_RefGene_Name's

hm450.reduced=hm450[-which(is.na(ucsc_ref)),]

dim(hm450)
dim(hm450.reduced)

#extract the names of the remaining probes after removing na entries of UCSC Ref Gene Name
Mred.colnames=colnames(M)[-which(is.na(ucsc_ref))]

#extract hm 450 UCSC ref gene names and parse by ";" into list
UCSC_ref_gn_parsed=strsplit(hm450.reduced$UCSC_RefGene_Name,";")

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
triobuildlist=match.gn.parsed(express.genenames, methyl.gn.parsed)

#get a final vector of all indicies pertaining to expression matched probes
Express.matched.methyl.probes=Mred.colnames[unqiue(unlist(triobuildlist))]



#check to make sure all rows agree:
all.equal(row.names(E), row.names(M))

#----------------------------------------------------------------------------------------------
#ASSEMBLING THE 4 DATA FILES (E, M, M_metadata, and E_metadata)

#match Expression data to Illumina ID's in biomart:
matching.E=match(colnames(E), gb37$IlmnID)

#make sure they match
print(head(gb37[na.omit(matching.E),]))
print(head(colnames(E)[-which(is.na(matching.E))]))

E.IlmnID.matching=E[-which(is.na(matching.E))]
E_metadata=gb37[na.omit(matching.E),]












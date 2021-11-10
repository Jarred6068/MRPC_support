


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
M=loadRData(fileName="/mnt/ceph/jarredk/Methyl/MethylData/MethylData.RegressResids.asMatrix.Rdata")
#read in the UCSC ref genenames for Methylation Data
M.genenames=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/MethylData/MethylData.RegressResids.UCSC_RefGene_Name.Rdata")
#Read in Human Methylation 450 data
hm450_meta=read.csv(file="/mnt/ceph/jarredk/Methyl/GPL13534_HumanMethylation450_15017482_v.1.1.csv")
#----------------------------------------------------------------------------------------------
#FINDING MATCHING METHYLATION AND EXPRESSION PROBES

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

save(UCSC_ref_gn_parsed, file = "/mnt/ceph/jarredk/Methyl/UCSC_ref_gn_parsed.Rdata")


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

#simple switch 
run=0
if(run==1){
  #a list of the methylation probes matched to each gene in express.genenames
  triobuildlist=match.gn.parsed(E.genenames, UCSC_ref_gn_parsed)

  save(triobuildlist, file = "/mnt/ceph/jarredk/Methyl/ExpressData/triobuildlist.Rdata")
}else{
  triobuildlist=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/ExpressData/triobuildlist.Rdata")
}

fn=function(x){x=ifelse(length(x)==0,NA,x)}

tbl.int0.removed=lapply(triobuildlist, fn)

#get a final vector of all indicies pertaining to expression matched probes
Express.matched.methyl.probes=Mred.colnames[unqiue(na.omit(unlist(tbl.int0.removed)))]

#save this list of probes to ExpressData folder
save(Express.matched.methyl.probes, file = "/mnt/ceph/jarredk/Methyl/ExpressData/Express.matched.methyl.probes.Rdata")

#----------------------------------------------------------------------------------------------
#matching up rows between methylation and expression

#check to make sure all rows agree:
rnM=NULL
for(i in 1:length(row.names(M))){rnM[i]=strsplit(row.names(M)[i], "/")[[1]][1]}
all.equal(row.names(E), rnM)

#match up rows and align expression to methylation
rows_matched=match(rnM, row.names(E))
print(sum(is.na(rows_matched)))
E.aligned=E[rows_matched,]

#save the aligned methylation data:
save(E.aligned, file = "/mnt/ceph/jarredk/Methyl/ExpressData/Expression_data_BSGSrows_aligned.Rdata")



#----------------------------------------------------------------------------------------------
#ASSEMBLING THE 4 DATA FILES (E, M, M_metadata, and E_metadata)


loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#read in BioMart Data containing Illumina ID's
gb37=read.table(file="/mnt/ceph/jarredk/Methyl/mart_export_GB37.txt", header = T, sep = ",")
#rename cols
colnames(gb37) = c("Gene.name", "chr", "start.bp", "end.bp", "IlmnID")

E.aligned=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/ExpressData/Expression_data_BSGSrows_aligned.Rdata")
E.genenames=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/ExpressData/Peerdata.bsgs_GeneIDs.Rdata")
#match Expression data to Illumina ID's in biomart:
#matching.E=match(colnames(E.aligned), gb37$IlmnID)
#get info for genes matching IlmnID and Gene name! (because there are occasionally diff genes sharing an ID)
get.match=function(E=NULL, E.gn=NULL, gb37=NULL){
  idx=NULL
  cn=colnames(E)
  for(i in 1:length(E.gn)){
    
    #if more than two genes sharing this ILMN ID then we include gene.name into search
    if(length(unique(gb37[which(gb37$IlmnID==cn[i]),]$Gene.name))>=2){
      
      ck=which(gb37$Gene.name==E.gn[i] & gb37$IlmnID==cn[i])
      #print(ck)
      
      if(length(ck>0)){
        idx[i]=ck[1]
      }else{
        idx[i]=NA
      }
      
    }else{
      
      ck=which(gb37$IlmnID==cn[i])
      #print(ck)
      
      if(length(ck>0)){
        idx[i]=ck[1]
      }else{
        idx[i]=NA
      }
      
    }
    
  }
  
  return(idx)
}

#matching.E=match(E.genenames, gb37$Gene.name)

matching.E=get.match(E=E.aligned, E.gn=E.genenames, gb37=gb37)

#make sure they match
print(head(gb37[na.omit(matching.E),]))
print(head(colnames(E.aligned)[-which(is.na(matching.E))]))

#construct the metadata file for expression after alignment by Illumina ID from BioMart
E.IlmnID.matching=E.aligned[,-which(is.na(matching.E))]
E.IlmnID.matching.gn=E.genenames[-which(is.na(matching.E))]
E_metadata=gb37[na.omit(matching.E),]

E_metadata=cbind.data.frame(ID=E_metadata$IlmnID, 
                            chr=E_metadata$chr, 
                            coordinate=E_metadata$start.bp, 
                            Gene.name=E_metadata$Gene.name)
A=c(1:22,"X","Y")
#remove rows with chr numbers not in 1:22, X, Y
idx2=which(E_metadata$chr %in% A)
E_metadata2=E_metadata[idx2,]
E.IlmnID.matching2=E.IlmnID.matching[,idx2]
E.IlmnID.matching.gn2=E.IlmnID.matching.gn[idx2]

#a function to check and make sure the E_metadata matches the Edata files before saving: if correct all TRUE
check.match=function(E=NULL, E.gn=NULL, E.meta=NULL){
  
  cn=colnames(E)
  TF.vec=NULL
  for(i in 1:dim(E_metadata)[1]){
    
    #print(isTRUE(cn[i]==E.meta$ID[i] & E.gn[i]==E.meta$Gene.name[i]))
    TF.vec[i]=ifelse(isTRUE(cn[i]==E.meta$ID[i] & E.gn[i]==E.meta$Gene.name[i]),i,NA)
  
  }
  return(TF.vec)
}


#some of the gene names won't match as the ones attached to bsgs data are outdated 
do.they.match=check.match(E=E.IlmnID.matching2, E.gn = E.IlmnID.matching.gn2, E.meta=E_metadata2)
head(E_metadata[which(is.na(do.they.match)),])
print(paste("Are All ILMN ID's Matching In Metadata?: ", all.equal(colnames(E.IlmnID.matching), E_metadata$ID), sep=""))


#save the metadata file for expression
save(E_metadata2, file = "/mnt/ceph/jarredk/Methyl/ExpressData/Express.BM.aligned.Ilmn.aligned.MetaData.Rdata")
#save the expression data containing only the probes matched to illumina entries
save(E.IlmnID.matching2, file = "/mnt/ceph/jarredk/Methyl/ExpressData/Express.BM.aligned.Ilmn.aligned.Rdata")
#save the attached bsgs gene names reduced
save(E.IlmnID.matching.gn2, file = "/mnt/ceph/jarredk/Methyl/ExpressData/Express.BM.aligned.Ilmn.aligned.Genenames.Rdata")

#now construct the metadata file for Methylation Probes
#since we already aligned hmm450 data to Methylation data we need only select the relevant columns
M_metadata=cbind.data.frame(IlmnID=hm450$IlmnID, chr=hm450$CHR, coordinate=hm450$MAPINFO)

save(M_metadata, file = "/mnt/ceph/jarredk/Methyl/ExpressData/All.Mprobes.Metadata.Rdata")









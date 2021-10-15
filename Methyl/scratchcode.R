







data1=cbind.data.frame(tmethyl.Mlt[,1], covariates.final)
colnames(data1)=c("log.Y", colnames(covariates.final))

model=lm(log.Y ~ Gender*poly(Age, 2), na.action = na.exclude, data = data1)








qc=function(X, V){ X==V }

match.gn.parsed=function(egn, mgnp){
  
  indexlist1=vector("list", length(egn))
  vec1=c(1:length(mgnp))
  timevec=NULL
  vec11=NULL
  
  for(i in 1:length(egn)){
    
    start.time=Sys.time()
    
    loglist1=lapply(mgnp, qc, egn[i])
    
    vec11=vec1[unlist(lapply(loglist1, any))]
    
    indexlist1[[i]]=vec11
    end.time=Sys.time()
    #print(vec11)
    print(end.time-start.time)
    timevec[i]=end.time-start.time
  }
  
  return(timevec)
  
}

est.time=match.gn.parsed(express.genenames[1:100], methyl.gn.parsed)
print(mean(est.time))







load(file="/mnt/ceph/jarredk/Methyl/Wksp1.Rdata")

fake_genos=function(csize=NULL, org_genos, expmat=NULL, expnames=NULL){
  percent.type=as.data.frame(matrix(0, nrow = dim(org_genos)[2], ncol = 4))
  
  for(i in 1:dim(org_genos)[2]){ 
    
    org_genos[,i]=as.factor(org_genos[,i]) 
    percent.type[i,]=summary(org_genos[,i])/sum(summary(org_genos[,i]))
    
  }
  
  avg.per=colMeans(percent.type)
  genos.mat=as.data.frame(matrix(0, nrow = csize, ncol = dim(expmat)[1] ))
  
  for(i in 1:csize){
    
    genos.mat[i, ]=sample(c(0,1,2,NA), dim(genos.mat)[2], replace = TRUE, prob = avg.per)
    
  }
  
  colnames(genos.mat)=row.names(expmat)
  row.names(genos.mat)=paste0("SNP", c(1:csize))
  return(t(genos.mat))
}

mat=fake_genos(1000000, ex.genotypes, express.aligned, express.genenames)
mat2=fake_genos(1000000, ex.genotypes, express.aligned, express.genenames)
mat3=fake_genos(1000000, ex.genotypes, express.aligned, express.genenames)
mat4=fake_genos(1000000, ex.genotypes, express.aligned, express.genenames)

genos.mat=cbind.data.frame(mat,mat2,mat3,mat4)

save(genos.mat, file = "/mnt/ceph/jarredk/Methyl/fakegenosBIG.Rdata")











#histograms on original data to check clustering



tmethyl.MO = loadRData(fileName = "/mnt/ceph/jarredk/Methyl/Mdata.matched.final.Rdata")

colstoplot=c("cg01516881","cg02798874","cg12213037","cg15383120","cg18110333","cg25004176","cg25817503","cg06758191","cg19601530","cg26668828")

matches=match(colstoplot, colnames(tmethyl.MO))

submat.methyl=tmethyl.MO[-1, matches ]

for(i in 1:dim(submat.methyl)[2]){
  
  png(paste("/mnt/ceph/jarredk/Methyl/ClusCheckPlots2/","Original_plots", colnames(submat.methyl)[i], ".png", sep = ""))
  hist(as.numeric(submat.methyl[,i]), 
       breaks = 20,
       xlab = "Methylation values",
       main = paste("Histogram:",colnames(submat.methyl)[i], sep = " "))
  
  dev.off()
  
}











#retrieve the probe ID's from the complete trio list obj

probenames=list()

cn=function(X){ v=colnames(X)[2]; return(v)}

samples=NULL

for(i in 1:length(EM.triolist)){
  
  samples[i]=length(lapply(EM.triolist[[i]], length))
  
}



for(i in 1:length(EM.triolist)){
  
  if(is.null(unlist(lapply(EM.triolist[[i]], cn)))){
    
    probenames[[i]]=NA
    
  }else{
    
    probenames[[i]]=unlist(lapply(EM.triolist[[i]], cn))
    
  }
  
}



genos.mat=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/fakegenosBIG.Rdata")
EM.triolist=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/completedTrios1.Rdata")
triobuildlist=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/triobuildlist.Rdata")
simulated.meta=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/fakegenoMeta.Rdata")
Mprobenames.final=loadRData(fileName = "/mnt/ceph/jarredk/Methyl/Used.Mprobes.final.Rdata")
meta_M=read.csv(file = "/mnt/ceph/megheib/M_G_data/GPL13534_M.csv")
#retreive relevant metadata for the "In-Use" methylation probes:
imp.meta_M=meta_M[,c(1,12,13)]
matched.meta_M=match(Mprobenames.final, imp.meta_M[,1])
imp.meta_M.final=imp.meta_M[matched.meta_M,]





#testing code for calc.cors
snps.in.range=list()
chr.GM.matched=imp.meta_M.final[which(imp.meta_M.final[,2]=="1"),]
chr.geno.matched=simulated.meta[which(simulated.meta[,1]=="1"),]
snps.in.range[[1]]=which(abs(chr.GM.matched[1,3]-chr.geno.matched[,2])<1000000)

cols=match(names(try1$Mprobes[[1]]), colnames(methyl.resids2))
methyl.mat=methyl.resids2[,cols]
cor_j=cor(methyl.mat[,1], genos.mat[, try1$Mprobes[[1]][[1]] ], use = "complete.obs")
names(cor_j)=colnames(genos.mat[,try1$Mprobes[[1]][[1]] ])



























total.value = sum(sum(G.cont)+sum(Busi.match)+sum(return.on.investment))
print(total.value)






























############################################
# R script for calculating correlations
# for all the genotypes
# 
# Note: Please change the file path and provide the filenames for 
# the genotype data and the genotype meta data where necessary.
#
# Note: all input files must be in the same folder specified in /filepath/Input
# Also, the user must create a folder called "Output" in the /filepath/ for the
# output tables
#
# The script currently runs each chromosome sequentially, from chr1 to chr22, plus chrX and chrY.
# If genotype data are not available for chrX and chrY, please remove those two lines (starting at line 306)
# 
# Alternatively, this script may be duplicated for each chromosome 
# such that the calculation may be performed in parallel.
#
#############################################


#####################
#set the loading and saving locations and bp proximity
#####################
# specify the size of the neighborhood
bp.set.range=500000
# specify the path to the input files and path to the output files
filepath = "/mnt/ceph/jarredk/Methyl/Correlation_Calculation/"
input.location = paste (filepath, "Input/", sep='')
# the user must create the folder "Output" within filepath/
output.location = paste (filepath, "Output/", sep='')
# This example generates two output files: filename1 and filename2
# The filenames are automatically generated by function calc.corsV2 below 

# load functions in calcCOR_v2.R for correlation calculation
source(file = paste(filepath, "calcCOR_v2.R", sep = ""))

####################
#load in Data
# Note that only the methylation data file (Mdata.Rdata) takes ~30 sec to load
# Other files are loaded quickly
# the user MUST fill in the filename for the genotype data and genotype metadata files
####################
# genotype data: a data matrix of 606 individuals in the rows and 1000 SNPs in the columns
genotype_data=loadRData(file = "/mnt/ceph/jarredk/Methyl/fakegenosBIG_2.Rdata")
# genotype metadata: a meta data matrix of 1000 SNPs in the rows and 2 columns [chr  coordinate]
G_metadata=loadRData(file = "/mnt/ceph/jarredk/Methyl/fakegenoMeta_2.Rdata")
# Methylation metadata: a meta data matrix of 363,516 probes and 3 columns [Ilmn_ID    chr    coordinate]
M_metadata=loadRData(fileName = paste(input.location, "meta_M.Rdata", sep = ""))
# Expression metadata: a meta data matrix of xxx probes and 3 columns [gene_ID     chr       coordinate]
E_metadata=loadRData(fileName = paste(input.location, "meta_E.Rdata", sep = ""))
# Expression Data: a matrix of 606 individuals in rows and xxx genes/probes in columns
Edata=loadRData(fileName = paste(input.location, "Expression_data_BM_aligned.Rdata", sep = ""))
# Methylation Data: a matrix of 606 individuals in rows and 363,516 probes in columns
Mdata=loadRData(fileName = paste(input.location, "Mdata.Rdata", sep = ""))


######################
#Run for each Chr
# Chr # is specified in the 'chrs' argument
######################

calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data,
            GMInfo = M_metadata,
            GEInfo = E_metadata,
            genoInfo = G_metadata,
            chrs = c("1"),
            fn = output.location,
            bp.range = bp.set.range)



#Run for Chr 2
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("2"),
            fn=output.location,
            bp.range = bp.set.range)


#Run for Chr 3
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("3"),
            fn=output.location,
            bp.range = bp.set.range)

#Run for Chr 4
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("4"),
            fn=output.location,
            bp.range = bp.set.range)








































#get the idx's for each mat for chr1 
Gidx=which(G_metadata[,2]=="1")
Midx=which(M_metadata[,2]=="1")
Eidx=which(E_metadata[,2]=="1")

#reduce meta data matrices
chr.geno.matched=G_metadata[Gidx,]
chr.GM.matched=M_metadata[Midx,]
chr.GE.matched=E_metadata[Eidx, ]




















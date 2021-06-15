


#A fast example using fake genotype data

#set the loading and saving locations and bp proximity

bp.set.range=500000
load.location="/mnt/ceph/jarredk/Methyl/Correlation_Code/"
save.location="/mnt/ceph/jarredk/Methyl/Correlation_Code/"

#load in search functions:
load(paste(load.location, "functions.Rdata", sep = ""))
#load in Data
genotype_data=read.table(file=paste(load.location, "EX_genotype_data.txt"), sep="\t", header=T)
G_metadata=read.table(file = paste(load.location, "EX_Genotype_Metadata.txt", sep = ""), sep = "\t", header = T)
M_metadata=loadRData(fileName = paste(load.location, "meta_M.final.Rdata", sep = ""))
E_metadata=loadRData(fileName = paste(load.location, "meta_E.final.Rdata", sep = ""))
Edata=loadRData(fileName = paste(load.location, "Expression_data_BM_aligned.final.Rdata", sep = ""))
Mdata=loadRData(fileName = paste(load.location, "Mdata.final.Rdata", sep = ""))


#Run for Chr 1
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("1"),
            fn=save.location,
            bp.range = bp.set.range)



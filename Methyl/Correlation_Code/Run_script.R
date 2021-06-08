
#set the loading and saving locations and bp proximity

bp.set.range=500000
load.location="/mnt/ceph/jarredk/Methyl/Correlation_Code/"
save.location="/mnt/ceph/jarredk/Methyl/cor_Lists_and_Tables/"

#load in search functions:
load(paste(load.location, "functions.Rdata", sep = ""))
#load in Data
genotype_data=loadRData(fileName = paste(load.location, "fakegenosBIG.Rdata", sep = ""))
G_metadata=loadRData(fileName = paste(load.location, "fakegenoMeta.Rdata", sep = ""))
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



#Run for Chr 2
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("2"),
            fn=save.location,
            bp.range = bp.set.range)


#Run for Chr 3
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("3"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 4
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("4"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 5
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("5"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 6
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("6"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 7
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("7"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 8
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("8"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 9
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("9"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 10
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("10"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 11
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("11"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 12
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("12"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 13
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("13"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 14
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("14"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 15
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("15"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 16
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("16"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 17
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("17"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 18
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("18"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 19
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("19"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 20
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("20"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 21
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("21"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr 22
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("22"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr X
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("X"),
            fn=save.location,
            bp.range = bp.set.range)

#Run for Chr Y
calc.corsV2(mmat = Mdata,
            emat = Edata,
            gmat = genotype_data, 
            GMInfo = M_metadata, 
            GEInfo = E_metadata, 
            genoInfo = G_metadata, 
            chrs=c("Y"),
            fn=save.location,
            bp.range = bp.set.range)



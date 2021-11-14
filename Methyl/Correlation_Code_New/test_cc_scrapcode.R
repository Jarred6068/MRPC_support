genotype_data2=cbind.data.frame(colnames(genotype_data)[-1], t(genotype_data[,-1]))
colnames(genotype_data2)=c("BGSG_ID", genotype_data[,1])


calc.cors(mmat = Mdata,
          emat = Edata,
          gmat = genotype_data2, 
          GMInfo = M_metadata, 
          GEInfo = E_metadata, 
          genoInfo = G_metadata, 
          chrs=c("1"),
          fn=output.location,
          bp.range = bp.set.range,
          snps.as.rows=FALSE)
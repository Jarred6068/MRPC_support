library(MRGN)

library(writexl)

data = loadRData('C:/Users/Bruin/OneDrive/Documents/GitHub/MRPC_support/HiC/sigfibro_tissue_interactions.RData')
master = loadRData('C:/Users/Bruin/OneDrive/Documents/GitHub/MRPC_support/Manuscript/Tables_extra/All_Trios_Master_table_LOND.RData')
# 
# 
# #write_xlsx(LOND, path = "/mnt/ceph/jarredk/Manuscript/LOND_HiC_master.xlsx")
# 
# tissnames = names(data)
# 
# 
# final.list = list()
# 
# for(i in 1:4){
#   
#   tabular = data[[i]]$summary.table1
#   final.list[[i]] = tabular
#   
# }
# 
# names(final.list) = tissnames
# 
# write_xlsx(final.list, path = "C:/Users/Bruin/OneDrive/Documents/GitHub/MRPC_support/HiC/LOND_Specificity_tables.xlsx")
# 
# 



create_table = function(data){
  tgid = data$CellsCulturedfibroblasts$summary.table1$trans.gene.ID[1:10]
  tgid.names = data$CellsCulturedfibroblasts$summary.table1$trans.gene.name[1:10]
  subtab = data$CellsCulturedfibroblasts$summary.table1[1:10, ]
  
  ct = NULL
  tisslist = NULL
  
  
  print(subtab)
  
  for(i in 1:10){
    
    tiss.ct = NULL
    print(paste0('Counting tissues for trans gene: ', tgid.names[i]))
    for(j in 1:4){
      print(paste0('tissue: ', names(data)[j]))
      ix = which(data[[j]]$summary.table1$trans.gene.ID == tgid[i])
      q = data[[j]]$summary.table1$qvals[ix]
      print(q)
      tiss.ct[j] = ifelse(q < 0.2, names(data)[j], NA)
      
    }  
    
    print(tiss.ct)
    
    tisslist[i] = paste(na.omit(tiss.ct), collapse = ',')
    ct[i] = length(na.omit(tiss.ct))
    
  }
  
  subtab$`# of Tissues Found` = ct
  subtab$`Tissues Found` = tisslist
  
  write.csv(subtab, file = 'C:/Users/Bruin/OneDrive/Documents/GitHub/MRPC_support/HiC/LOND_top_10_sig.csv')
}

#debug(create_table)
create_table(data)








































# Extract TF-target pairs from a dataset of trios
#
# input
# trios: table of trios.  Results from MRGN
# tf.name: a vector of TF names
# is.cis: logical.  TRUE if the cis gene is the TF
extractTFsFromTrios <- function (trios, tf.name, is.cis = TRUE) {
  tfs.gtex <- NULL
  for (i in 1:length(tf.name)) {
    print (i)
    
    # look for TF in trios
    tmp.match <- ifelse(is.cis, which (trios$Cis.Gene.Name==tf.name[i]), which (trios$Trans.Gene.Name==tf.name[i]))
    
    # if TF is found
    # extract all the rows with different targets
    if (length(tmp.match) & !is.na(tmp.match)) {
      for (j in 1:length (tmp.match)) {
        tfs.gtex <- rbind (tfs.gtex, c(trios[tmp.match[j],c(1, 3, 4,9, 10,23:25)]))
        #print (trios.all.m12[which (trios.all.m12$Trans.Gene.Name==tf.name[i]),c(3,9,22:24)])
      }
    } 
  }
  
  return (tfs.gtex)
}


# function to extract TF-target pairs from all tissues in GTEx
# input
# tf.name: a vector of TF names
# tf.targets: a list.  Each element is a vector of target gene names
# tissues: a vector of tissue names
# file.path: location of GTEx RData objects
# db.entrez: list of Entrez IDs for gene names
# db.ensembl: list of ENSEMBL IDs for Entrez IDs
#
# output
# data table.  Each row is a TF-target pair found in GTEx trio analysis.
extractTFTargetFromGTEx <- function (tf.name, tf.targets, tissues, file.path, db.entrez, db.ensembl) {
  result <- NULL
  # read in each tissue in a for loop
  for (i in 1:nrow(tissues.vec)) {
    #for (i in 1:5) {
    print (i)
    file1 <- paste(file.path,
                   tissues[i],"_AllPC/data.snp.cis.trans.final.",
                   tissues[i],".V8.unique.snps.RData", sep = "")
    data <- loadRData(fileName=file1)
    
    print (tissues[i])
    
    # locate columns in GTEx using the TF Ensembl IDs
    # look for targets around the TF using target Ensembl IDs
    # columns of the TF +/-1 should help locate the gene in the same trio
    # if trios are found, compute the cor between the two genes
    # return the tissue, two genes and their correlation
    for (k in 1:length (tf.name)) {
      #  for (k in 1:5) {
      cat (k)
      # extract Entrez IDs for each TF-target set
      gene.entrez <- db.entrez[match(c(tf.name[k], tf.targets[[k]]), names (db.entrez))]
      # extract Ensembl IDs using Entrez IDs
      gene.ensembl <- db.ensembl[unlist (gene.entrez)]
      
      tf.ensembl.gtex <- na.omit (match (unlist (gene.ensembl)[1], colnames(data)))
      if (length (tf.ensembl.gtex)>0) {
        # extract the neighbors of the TF
        tmp.cols <- colnames (data)[c(tf.ensembl.gtex-1, tf.ensembl.gtex+1)]
        # look for targets among the neighbors
        targets.ensembl.gtex <- tmp.cols[na.omit (match (unlist(gene.ensembl)[-1], tmp.cols))]
        # if found, 
        if (length(targets.ensembl.gtex)) {
          for (j in 1:length (targets.ensembl.gtex)) {
            # 
            tmp.cor <- cor (data[,tf.ensembl.gtex[1]], data[,match (targets.ensembl.gtex[j], colnames(data))])
            result <- rbind (result, c(tissues[i], unlist(gene.ensembl)[1], tf.name[k], targets.ensembl.gtex[j], tmp.cor))
            print (result)
          }
        }
      }
    }
    print("\n")
  }
  
  return (result)
}
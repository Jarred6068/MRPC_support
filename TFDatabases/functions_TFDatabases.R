#######################################
# function to load RData and rename it
#######################################
loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


###############################################################################
# function to read in multiple data tables of trios, each table from a tissue
###############################################################################
readMRGNTrios <- function (filename) {
  sheets <- openxlsx::getSheetNames(filename)
  trios <- lapply(sheets, openxlsx::read.xlsx, xlsxFile=filename)
  # assigning names to data frame
  names(trios) <- sheets
  ntissues <- length (trios)
  result <- rbindlist (trios, use.names = TRUE, fill = TRUE, idcol = "Tissue")
  
  return (result)
}


####################################################
# Extract TF-target pairs from a dataset of trios
#
# input
# trios: table of trios.  Results from MRGN
# tf.name: a vector of TF names
# is.cis: logical.  TRUE if the cis gene is the TF
####################################################
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



#######################################################################
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
# data table: Each row is a TF-target pair found in GTEx trio analysis.
#######################################################################
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
        # if found, compute the cor between the two genes
        # return the tissue, two genes and their correlation
        if (length(targets.ensembl.gtex)) {
          for (j in 1:length (targets.ensembl.gtex)) {
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


#######################################################################
# function to compare two tables
# Find same TF-target pairs in GTEx trios and the TF database
#######################################################################
compareTables <- function (trios.final, trios.candidate, is.cis=TRUE) {
  if (is.cis) {
    # if cis gene is the mediator
    # find TF-target pairs in both sets
    #    for (i in 1:nrow(trios.candidate)) {
    #      print (i)
    #      tmp <- which (trios.final[,1]==trios.candidate[i,2] & trios.final[,3] == trios.candidate[i,4])
    #      print (trios.final[tmp,])
    #    }
    
    for (i in 1:nrow(trios.final)) {
      print (i)
      tmp <- which (trios.candidate[,2] == trios.final[i,2] & trios.candidate[,4] == trios.final[i,4])
      if (length(tmp)) {
        #        print (trios.candidate[tmp,])
        print (tmp)
        print (unlist (trios.final[i,]))
      }
    }
  } else {
    #    for (i in 1:nrow(trios.candidate)) {
    #      print (i)
    #      tmp <- which (trios.final[,3]==trios.candidate[i,2] & trios.final[,1] == trios.candidate[i,4])
    #      print (trios.final[tmp,])
    #    }
    
    for (i in 1:nrow(trios.final)) {
      print (i)
      tmp <- which (trios.candidate[,2] == trios.final[i,4] & trios.candidate[,4] == trios.final[i,2])
      if (length(tmp)) {
        #       print (trios.candidate[tmp,])
        print (tmp)
        print (unlist (trios.final[i,]))
      }
    }
    
  }
  
}


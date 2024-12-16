# function to compare two tables
# Find same TF-target pairs in GTEx trios and the TF database
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

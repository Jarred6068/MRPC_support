# function to read in multiple data tables of trios, each table from a tissue
readMRGNTrios <- function (filename) {
  sheets <- openxlsx::getSheetNames(filename)
  trios <- lapply(sheets, openxlsx::read.xlsx, xlsxFile=filename)
  # assigning names to data frame
  names(trios) <- sheets
  ntissues <- length (trios)
  result <- rbindlist (trios, use.names = TRUE, fill = TRUE, idcol = "Tissue")
  
  return (result)
}

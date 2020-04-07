# UNDER CONSTRUCTION
# mapping algorithm for mapping a taxonomy table onto a taxonomy by exact name-matching, regardless of rank

# taxin is a dataframe containing ASV sequences and a set of taxonomic assignments
# tax2map2 can either be a filename or a dataframe of the taxonomic nomenclature you'd like to "translate" taxin into (use pr2_all_tax.csv)
# synonym.file is a filename or dataframe (use tax_synonyms_FINAL.csv in this folder to build) containing synonyms to check if exact name matching doesn't work
# outfilez is a character vector specifying names of output csv files. if none, do not write outputs to files

# function should always return a dataframe and a character vector, and potentially a second dataframe
# the first dataframe should include all unique taxonomic assignments (unique rows of taxin, excluding ASV seqs) and their corresponding mapped assignments
# the character vector should include all names within taxin that were unable to be mapped
# the optional second data frame should have the ASV seqs supplied in taxin with the mapped taxonomic assignments (structure should be identical to taxin except for taxonomic nomenclature)

taxmapper <- function(taxin, tax2map2, 
                      synonym.file = "tax_synonyms_FINAL.csv", 
                      outfilez = "none") {
  ref_cols <- rev(colnames(taxin))
  map_cols <- rev(colnames(tax2map2))
  
  ref_mapped <- NA
  mapped <- FALSE
  ref_not_mapped <- vector()
  
  for (row in 1:nrow(taxin)) {
    for (col in 1:ncol(taxin)) {
      if (!is.null(taxin[row, col])) {
        result <- findRow(row[, ref_cols[col]], tax2map2, map_cols)
        if (!is.empty(result) & length(result) == 1) {
          ref_mapped[row] <- result
          mapped <- TRUE
          break
        }
        else {
          if (row[colNames(taxin)[1]] != 'Bacteria' & row[colNames(taxin)[1]] != 'Archaea' & !mapped) {
            ref_not_mapped <- c(ref_not_mapped, taxin[row, ref_cols[col]])
          }
        }
      }
      if (row[colNames(taxin)[1]] != 'Bacteria' & row[colNames(taxin)[1]] != 'Archaea' & !mapped) {
        ref_mapped[row] <- data.frame(c('Bacteria', rep(NA, times = length(tax2map2))), colNames = colNames(tax2map2))
        mapped <- TRUE
      }
      if (!mapped) {
        ref_not_mapped <- c(ref_not_mapped, row)
      }
      else {
        mapped <- FALSE
      }
    }
  }
  
  ref_mapped_index <- vector(ref_mapped.keys)
  
  ref_mapped_df <- data.frame(colnames = colNames(taxin))
  map_mapped_df <- data.frame(colnames = colNames(tax2map2))
  
  for (index in ref_mapped_index) {
    ref_mapped_df <- rbind(ref_mapped_df, ref_df[index,])
    map_mapped_df <- rbind(map_mapped_df, ref_mapped[index])
  }
  
  combined <- rbind(ref_mapped_df, map_mapped_df)
  
  return (combined)
}




findRow <- function(colName, dataFrame, listOfCols) {
  for (col in 1:length(listOfCols)) {
    result <- dataFrame[, dataFrame$listOfCols[col] == colName]
    if (length(result) > 0) {
      mapped_row <- rep(NA, times = length(listOfCols))
      mapped_row[1:length(listOfCols)-col] = result[0,1:length(listOfCols)-col]
      return (data.frame(mapped_row, colName = rev(listOfCols)))
    }
    return (data.frame())
  }
}
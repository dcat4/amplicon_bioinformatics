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

# helper function takes in given taxonomy and finds a potential matching row from tax2map2
# returns the matching row from tax2map2 or NA if it can't find any

findMapping <- function(taxonomy, tax2map2) {
  cols <- rev(names(tax2map2))
  for (i in 1:length(cols)) {
    matchings <- tax2map2[which(tax2map2[, cols[i]] == taxonomy), ]
    if (nrow(matchings) != 0) {
      matched.row <- data.frame(matrix(rep(NA, length(cols)), ncol = length(cols), nrow = 1))
      colnames(matched.row) <- rev(cols)
      matched.row[1:(length(cols)-i+1)] <- matchings[1, ][1:(length(cols)-i+1)]
      return (matched.row)
    }
  }
  return (NA)
}

taxmapper <- function(taxin, tax2map2, 
                      synonym.file = "tax_synonyms_FINAL.csv", 
                      outfilez = "none") {
  
  # remove duplicates and remove sVN and ASV columns
  taxin.u <- unique(taxin[,-c(1,2)])
  tax2map2.u <- unique(tax2map2[,-c(1,2)])
  
  taxin.cols <- rev(names(taxin.u))
  
  # non eukaryotes
  nonexist <- c('Bacteria', 'Archaea')
  
  # not mapped
  not.mapped <- vector()
  
  mapped <- data.frame(matrix(ncol=(ncol(taxin.u) + ncol(tax2map2.u)),nrow=0, dimnames=list(NULL, c(names(taxin.u), names(tax2map2.u)))))
  
  for (row in 1:nrow(taxin.u)) {
    for (col in 1:ncol(taxin.u)) {
      taxonomy <- taxin.u[row, taxin.cols[col]]
      if (!is.na(taxonomy)) {
        match <- findMapping(taxonomy, tax2map2.u)
        if (is.data.frame(match)) {
          combined <- cbind(taxin.u[row, ], match)
          mapped <- rbind(mapped, combined)
          break
        }
        else {
          if (is.element(taxonomy, nonexist)) {
            null.row <- data.frame(matrix(rep(NA, ncol(tax2map2)), ncol = ncol(tax2map2), nrow = 1, dimnames=list(NULL, names(tax2map2.u))))
            null.row[1] <- 'Bacteria'
            combined <- cbind(taxin.u[row, ], null.row)
            mapped <- rbind(mapped, combined)
          }
          else {
            not.mapped <- c(not.mapped, taxonomy)
          }
        }
      }
    }
  }
}


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

taxmapper <- function(taxin, tax2map2, exceptions,
                      synonym.file = "tax_synonyms_FINAL.csv", 
                      outfilez = "none") {
  
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
  
  getSynonyms <- function(taxonomy, syn.df) {
    found.rows <- syn.df[which(syn.df == taxonomy, arr.ind=TRUE)[,'row'],]
    if (length(found.rows) > 0) {
      v <- as.character(as.matrix(found.rows))
      return (unique(c(taxonomy, v[!is.na(v)])))
    }
    else {
      return (c(taxonomy))
    }  
  }
  
  taxin.u <- unique(taxin[, !(names(taxin) %in% c("svN","ASV"))])
  tax2map2.u <- unique(tax2map2)
  
  synonyms <- read.csv(synonym.file)
  synonyms <- synonyms[, colnames(synonyms)[startsWith(colnames(synonyms), "Name")]]
  
  taxin.cols <- rev(names(taxin.u))
  not.mapped <- vector()
  mapped <- data.frame(matrix(ncol=(ncol(taxin.u) + ncol(tax2map2.u)),nrow=0, dimnames=list(NULL, c(names(taxin.u), names(tax2map2.u)))))
  
  for (row in 1:nrow(taxin.u)) {
    highest.tax <- taxin.u[row, taxin.cols[ncol(taxin.u)]]
    if (is.element(highest.tax, exceptions)) {
      null.row <- data.frame(matrix(rep(NA, ncol(tax2map2.u)), ncol = ncol(tax2map2.u), nrow = 1, dimnames=list(NULL, names(tax2map2.u))))
      null.row[1] <- highest.tax
      combined <- cbind(taxin.u[row, ], null.row)
      mapped <- rbind(mapped, combined)
    }
    else {
      for (col in 1:ncol(taxin.u)) {
        tax <- taxin.u[row, taxin.cols[col]]
        pos.taxs <- getSynonyms(tax, synonyms)
        matched <- FALSE
        for (taxonomy in pos.taxs) {
          last <- FALSE
          if (match(taxonomy, pos.taxs) == length(pos.taxs)) {
            last <- TRUE
          }
          if (!is.na(taxonomy)) {
            match <- findMapping(taxonomy, tax2map2.u)
            if (is.data.frame(match)) {
              combined <- cbind(taxin.u[row, ], match)
              mapped <- rbind(mapped, combined)
              matched <- TRUE
              break
            }
            else {
              if (last) {
                not.mapped <- c(not.mapped, tax)
              }
            }
          }
        }
        if (matched) {
          break
        }
      }
    }
  }
  
  names(taxin) <- c(colnames(taxin)[1:2], toupper(colnames(taxin.u)))
  mapped.t <- mapped
  names(mapped.t) <- c(toupper(colnames(mapped)[1:ncol(taxin.u)]), colnames(mapped)[(ncol(taxin.u)+1):ncol(mapped)])
  
  asv.mapped <- merge(x=taxin, y=mapped.t, by=toupper(colnames(taxin.u)), all.x=TRUE)
  asv.mapped <- asv.mapped[ , !(colnames(asv.mapped) %in% toupper(colnames(taxin.u)))]
  
  not.mapped <- unique(not.mapped)
  
  if (outfilez != "none") {
    write.csv(mapped, outfilez[1])
    write.csv(as.data.frame(not.mapped), outfilez[2])
    write.csv(asv.mapped, outfilez[3])
  }
  
  return (list(mapped, not.mapped, asv.mapped))
}


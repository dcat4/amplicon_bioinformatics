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

taxmapper <- function(taxin, tax2map2, exceptions, ignore.format = FALSE,
                      synonym.file = "tax_synonyms_FINAL.csv", 
                      outfilez = "none") {
  
  # function to create alternative terms with the suffix
  createAlts <- function(taxs) {
    result <- vector()
    for (tax in taxs) {
      a1 <- gsub("(phyta)", "phyceae", tax)
      a2 <- gsub("(phyta)", "phyte", tax)
      a3 <- gsub("(phyta)", "phytes", tax)
      a4 <- gsub("(phyceae)", "phyta", tax)
      a5 <- gsub("(phyceae)", "phyte", tax)
      a6 <- gsub("(phyceae)", "phytes", tax)
      a7 <- gsub("(phyte)", "phyta", tax)
      a8 <- gsub("(phyte)", "phyceae", tax)
      a9 <- gsub("(phyte)", "phytes", tax)
      a10 <- gsub("(phytes)", "phyta", tax)
      a11 <- gsub("(phytes)", "phyceae", tax)
      a12 <- gsub("(phytes)", "phyte", tax)
      result <- c(result, unique(c(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12)))
    }
    return(unique(result))
  }
  
  # function to remove hyphens, underscores, upper case of name
  preprocessTax <- function(taxonomy) {
    # split terms by hyphens
    no.hyphen <- strsplit(taxonomy, "-")
    # split terms by underscores
    no.underscore <- strsplit(taxonomy, "_")
    # split terms by first instance of underscores and combine previous splits
    taxs <- c(no.hyphen[[1]], no.underscore[[1]], gsub("(_.*)", "", taxonomy), taxonomy)
    # remove duplicates
    taxs <- unique(taxs)
    # convert all to lowercase
    no.upper <- tolower(taxs)
    # create alternative suffixes for certain taxonomies
    final.taxs <- createAlts(unique(c(taxs, no.upper)))
    return(final.taxs)
  }
  
  # function to search through the tax2map2 to find a match for the taxonomy name inputted
  findMapping <- function(taxonomy, tax2map2) {
    # iterate through the most specific ranking to the most generic ranking
    cols <- rev(names(tax2map2))
    for (i in 1:length(cols)) {
      matchings <- tax2map2[which(tax2map2[, cols[i]] == taxonomy), ] # find rows that match at that rank
      if (nrow(matchings) != 0) {
        # create respective row for the match
        # make everything downstream of the rank found to be NA's
        matched.row <- data.frame(matrix(rep(NA, length(cols)), ncol = length(cols), nrow = 1))
        colnames(matched.row) <- rev(cols)
        # grab only the first match found
        matched.row[1:(length(cols)-i+1)] <- matchings[1, ][1:(length(cols)-i+1)]
        return (matched.row)
      }
    }
    return (NA)
  }
  
  # function to search through the synonyms data frame to find synonyms for given taxonomy name
  getSynonyms <- function(taxonomy, syn.df) {
    found.rows <- syn.df[which(syn.df == taxonomy, arr.ind=TRUE)[,'row'],] # find rows for synonym
    if (length(found.rows) > 0) {
      # populate the taxonomy with its synonyms
      v <- as.character(as.matrix(found.rows))
      return (unique(c(taxonomy, v[!is.na(v)])))
    }
    else {
      # if no synonyms found, just return the taxonomy
      return (c(taxonomy))
    }  
  }
  
  # rename tax2map2 columns for uniqueness
  colnames(tax2map2) <- paste("tax2map2", colnames(tax2map2), sep="_")
  
  # grab only the taxonomies part of the dataframes
  taxin.u <- unique(taxin[, !(names(taxin) %in% c("svN","ASV"))])
  tax2map2.u <- unique(tax2map2)
  
  # read in the synonyms file
  synonyms <- read.csv(synonym.file)
  synonyms <- synonyms[, colnames(synonyms)[startsWith(colnames(synonyms), "Name")]]
  
  taxin.cols <- rev(names(taxin.u))
  
  # keep track of the taxonomy names that are not mapped
  not.mapped <- vector()
  
  # finialized mapping table from taxin to tax2map2 with only the taxonomy names
  mapped <- data.frame(matrix(ncol=(ncol(taxin.u) + ncol(tax2map2.u)),nrow=0, dimnames=list(NULL, c(names(taxin.u), names(tax2map2.u)))))
  
  # iterate through each row and column of taxin data frame
  for (row in 1:nrow(taxin.u)) {
    # keep track of the most generic taxonomy name
    highest.tax <- taxin.u[row, taxin.cols[ncol(taxin.u)]]
    # see if it is in the exceptions to skip the row
    if (is.element(highest.tax, exceptions)) {
      # create a NA row assignment since part of exceptions
      null.row <- data.frame(matrix(rep(NA, ncol(tax2map2.u)), ncol = ncol(tax2map2.u), nrow = 1, dimnames=list(NULL, names(tax2map2.u))))
      null.row[1] <- highest.tax
      combined <- cbind(taxin.u[row, ], null.row)
      mapped <- rbind(mapped, combined)
    }
    else {
      for (col in 1:ncol(taxin.u)) {
        # keep track of the original taxonomy name 
        orig.tax <- taxin.u[row, taxin.cols[col]]
        # process the name to get alternatives by igorning its format
        if (ignore.format) {
          pos.taxs <- preprocessTax(orig.tax)
          for(tax in pos.taxs) {
            pos.taxs <- c(pos.taxs, getSynonyms(tax, synonyms))
          }
          pos.taxs <- unique(c(orig.tax, pos.taxs))
        }
        else {
          pos.taxs <- getSynonyms(orig.tax, synonyms)
        }
        # flag to keep track of when the row is already matched
        matched <- FALSE
        # counter to keep track what column number we are on
        counter <- 1
        # iterate through all alternatives of the taxonomy name with original one first
        for (taxonomy in pos.taxs) {
          last <- FALSE
          if (counter == length(pos.taxs)) {
            last <- TRUE
          }
          if (!is.na(taxonomy)) {
            # find matching
            match <- findMapping(taxonomy, tax2map2.u)
            if (is.data.frame(match)) {
              combined <- cbind(taxin.u[row, ], match)
              mapped <- rbind(mapped, combined)
              matched <- TRUE
              break
            }
            else { # if no matching is found, add to not.mapped
              if (last) {
                not.mapped <- c(not.mapped, orig.tax)
              }
            }
          }
          counter <- counter + 1
        }
        # if a match is found, we can move onto the next row
        if (matched) {
          break
        }
      }
    }
  }
  
  # use the mapped table created to left join the original data frame with ASV and svN
  asv.mapped <- merge(x=taxin, y=mapped, by=colnames(taxin.u), all.x=TRUE)
  asv.mapped <- asv.mapped[ , !(colnames(asv.mapped) %in% colnames(taxin.u))]
  # remove the unique addition of tax2map2 column names
  colnames(asv.mapped) <- gsub("tax2map2_", "", colnames(asv.mapped))
  
  # filter out duplicates for not mapped taxonomy names
  not.mapped <- unique(not.mapped)
  
  colnames(mapped) <- gsub("tax2map2_", "", colnames(mapped))
  
  if (outfilez != "none") {
    write.csv(mapped, outfilez[1], row.names=FALSE)
    not.mapped.df <- as.data.frame(not.mapped)
    write.table(not.mapped.df, outfilez[2], row.names=FALSE, col.names=FALSE)
    write.csv(asv.mapped, outfilez[3], row.names=FALSE)
  }
  
  return (list(mapped, not.mapped, asv.mapped))
}


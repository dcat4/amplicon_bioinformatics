# to improve the time complexity of the algorithm, try to map the unique values of taxin
# create a mapper dataframe so it can be left joined to the original avoiding additional searching
traitmapper_Ramond_fast <- function(taxin, map2,
                                    map2.taxnames = c("Lineage1","Lineage2","Lineage3","Lineage4",
                                                      "Lineage5","Lineage6","Lineage7","Lineage8",
                                                      "Lineage9","Lineage10","Lineage11","Lineage12",
                                                      "Fam","Taxogroup","Taxo1","Last"),
                                    dont.map = c("Eukaryota", "Archaea", "Bacteria", 
                                                 "Alveolata","Opisthokonta","Archaeplastida","Excavata","Rhizaria","Stramenopiles",
                                                 "Hacrobia","Amoebozoa","Apusozoa","Eukaryota_X","Protalveolata","Terrabacteria"),
                                    filezout = "none") {
  
  # finds the corresponding row in map2 to map to based on map2.taxnames
  # the row mapped will be the row with the least number of NA's
  # can be further sped up by arbitrarily taking the 1st row of the matchings
  findMapping <- function(taxonomy, map2, map2.taxnames) {
    cols <- rev(map2.taxnames)
    for (col in cols) {
      matchings <- map2[which(map2[, col] == taxonomy), ]
      if (nrow(matchings) != 0) {
        matched.row <- data.frame(matrix(NA, ncol = ncol(map2) + 2), nrow = 1)
        colnames(matched.row) <- c("name.mapped", "n.hitz", colnames(map2))
        matched.row$name.mapped <- taxonomy
        matched.row$n.hitz <- nrow(matchings)
        #matched.row[3:(ncol(map2) + 2)] <- matchings[which.max(rowSums(!is.na(matchings))), ]
        matched.row[3:(ncol(map2) + 2)] <- matchings[1,]
        return (matched.row)
      }
    }
    return (NA)
  }
  
  # get the unique taxonomies to map to
  taxin.u <- unique(taxin[, !(names(taxin) %in% c("svN", "ASV"))])
  taxin.cols <- rev(names(taxin.u)) # set up mapping by most specific grouping first
  not.mapped <- vector() # keep track of the taxonomies in taxin not able to map to
  
  # mapping data frame to fill in while iterating through the rows of unique taxin
  mapped <- data.frame(matrix(ncol=(ncol(taxin.u) + 2 + ncol(map2)), nrow=0, dimnames=list(NULL, c(names(taxin.u), "name.mapped", "n.hitz", names(map2)))))
  
  # create a NA row beforehand incase no matchings are found for a given row
  nan.row <- data.frame(matrix(NA, ncol = ncol(map2) + 2), nrow = 1)
  colnames(nan.row) <- c("name.mapped", "n.hitz", colnames(map2))
  
  for (row in 1:nrow(taxin.u)) {
    counter <- 1
    for (col in taxin.cols) {
      tax <- taxin.u[row, col]
      if (!is.na(taxonomy)) {
        match <- findMapping(tax, map2, map2.taxnames)
        if (is.data.frame(match)) {
          combined <- cbind(taxin.u[row, ], match)
          mapped <- rbind(mapped, combined)
          matched <- TRUE
          break
        } else {
          not.mapped <- c(not.mapped, tax)
        }
      }
      # no matchings found and finished searching most general column
      if (counter == length(taxin.cols)) {
        combined <- cbind(taxin.u[row, ], nan.row)
        mapped <- rbind(mapped, combined)
      }
      counter <- counter + 1
    }
  }
  
  # left join to the original taxin data frame with the mapped data frame
  taxin.mapped <- merge(x=taxin, y=mapped, by=colnames(taxin.u), all.x=TRUE)
  
  # put the svN and ASV cols in the front
  move <- c("svN", "ASV")
  others <- setdiff(names(taxin.mapped), move)
  others <- others[!is.na(want)] # remove unnecessary NA column for row nums
  taxin.mapped <- taxin.mapped[c(move, want)] # correct order of columns
  
  if (filezout != "none") {
    write.csv(taxin.mapped, filezout[1], row.names=FALSE)
    not.mapped.df <- as.data.frame(not.mapped)
    write.table(not.mapped.df, filezout[2], row.names=FALSE, col.names=FALSE)
  }
  
  return(list(taxin.mapped, not.mapped))
}



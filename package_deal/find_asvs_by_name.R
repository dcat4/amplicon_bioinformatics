# seems to work :)) needs more QC before generalizing

# trying to write something convenient to find ASVs that have a particular name across multiple tax tables
# that's also agnostic to taxonomic heirarchy (different numbers of ranks can be used in input tables)

# all tax tables should have the same ASVs as rows in the same order

# returns a dataframe where each column is one of names2find and contains row indices where the name
# was found across all input tax tables

# optionally will return a list with first element as above, and 
# 2nd element a list with each element a dataframe similar to first output but for each input tax table
# this list is in order of input tax tables

# ... = arbitrary number of taxonomy tables
# names2find = character vector of the names you want to find in each input tax table
# return.byTable = if TRUE, return a list with 2nd element a list (see above)

find_asvs_by_name <- function(..., names2find, return.byTable = FALSE) {
  
  xx <- list(...)
  
  for (i in 1:length(names2find)) {
    nn <- names2find[i]
    eh <- lapply(xx, FUN = function(x) apply(x, MARGIN = c(1,2), function(y) nn == y))
    matchr.list <- lapply(eh, FUN = function(x) which(rowSums(x, na.rm = TRUE) > 0))
  }
  
  # initialize list for storing index vectors by name for each tax table input:
  allout.list <- rep(list(base::as.data.frame(matrix(NA, nrow = nrow(xx[[1]]), ncol = length(names2find), dimnames = list(NULL,names2find)), stringsAsFactors = FALSE)), times = length(xx))
  # initialize dataframe for overlapping indices of each name across the input tax tables
  in.all <- as.data.frame(matrix(NA, nrow = nrow(xx[[1]]), ncol = length(names2find), dimnames = list(NULL,names2find)), stringsAsFactors = FALSE)
  for (i in 1:length(names2find)) {
    nn <- names2find[i]
    eh <- lapply(xx, FUN = function(x) apply(x, MARGIN = c(1,2), function(y) nn == y))
    matchr.list <- lapply(eh, FUN = function(x) which(rowSums(x, na.rm = TRUE) > 0))
    compme <- matchr.list[[1]] # for comparing to other list elements in matchr.list
    for (j in 1:length(matchr.list)) {
      plombus <- matchr.list[[j]]
      allout.list[[j]][1:length(plombus),nn] <- plombus # assign indices for this name-taxtable combo
      if (j > 1){
        # intersect the index vector to the first one
        # iteratively intersecting will trim it down
        compme <- intersect(compme, plombus)
      }
    }
    in.all[1:length(compme), nn] <- compme
  } 
  
  # remove rows of all NA:
  in.all <- in.all[!rowSums(is.na(in.all)) == ncol(in.all),]
  allout.list <- lapply(allout.list, FUN = function(x) subset(x, !rowSums(is.na(x)) == ncol(x)))
  
  if (return.byTable) {
    return(list(in.all, allout.list))
  } else {
    return(in.all)
  }
}
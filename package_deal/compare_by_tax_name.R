# inputs:
# 1. ... = arbitrary # of tax tables (may have to make this set)
# 2. tablenames = names of each tax table input
# 3. taxnames = the names for which you'd like to do the comparisons
# 4. return.conflix = TRUE or FALSE --> if TRUE, computes conflict and non-conflict indices and counts and reports them rather than detailed comparison options

# outputs:
# 1. a list of lists -- names of the list correspond to taxnames
# within each list the first element is a list including a summary dataframe (comparison options with counts)
# the second list is a dataframe of indices for each of the comparison options in the summary (columns corresponding to rows of summary)
# in other words, this statement should be TRUE: length(out[[2]][,1]) - length(which(is.na(out[[2]][,1]))) == out[[1]][1,"Freq"]
# 2. if return.conflix, the comparisons are condensed to 2 categories and return values are the same w/ only the 2 categories reported

compare_by_tax_name <- function(..., taxnames, tablenames, return.conflix = FALSE) {
  xx <- list(...)
  notuz <- nrow(xx[[1]])
  
  # initialize list for storing index vectors by name for each tax table input:
  allout.list <- vector(mode = "list", length = length(taxnames))
  # initialize dataframe for overlapping indices of each name across the input tax tables
  in.all <- as.data.frame(matrix(NA, nrow = nrow(xx[[1]]), ncol = length(taxnames), dimnames = list(NULL,taxnames)), stringsAsFactors = FALSE)
  for (i in 1:length(taxnames)) {
    nn <- taxnames[i]
    eh <- lapply(xx, FUN = function(x) apply(x, MARGIN = c(1,2), function(y) nn == y))
    # in here you need to determine the column of each tax table where the name is found, and compile other options for this rank
    colz <- lapply(eh, FUN = function(x) which(apply(x, MARGIN = 2, FUN = any, na.rm = TRUE)))
    look.here <- vector("list", length = length(xx))
    # return(list(look.here, colz, xx))
    for (j in 1:length(colz)){
      look.here[[j]] <- xx[[j]][,colz[[j]]]
    }
    # look.here includes the column of each taxtable where nn was found
    compme <- data.frame(matrix(unlist(look.here), nrow = notuz, ncol = length(xx), byrow = FALSE), stringsAsFactors = FALSE)
    colnames(compme) <- tablenames
    matchr.list <- lapply(eh, FUN = function(x) which(rowSums(x, na.rm = TRUE) > 0))
    hitterz <- unique(unlist(matchr.list)) # all rows where at least one taxonomy table contained the name
    ohyea <- as.data.frame(table(compme[hitterz,], useNA = "always"))
    ohyea <- ohyea[ohyea$Freq != 0,]
    ohyea.i <- as.data.frame(matrix(NA, nrow = notuz, ncol = nrow(ohyea))) # create an all NA dataframe of indices for each row of ohyea
    # change NAs to "NA" (character) so you they'll be ==
    compme[is.na(compme)] <- "NA"
    for (j in 1:(ncol(ohyea)-1)) {
      ohyea[,j] <- as.character(ohyea[,j])
    }
    ohyea[is.na(ohyea)] <- "NA"
    for (j in 1:nrow(ohyea)){
      this1 <- ohyea[j,!colnames(ohyea) %in% c("Freq")]
      this1.i <- which(apply(compme, 1, function(x) all(x == this1)))
      if (length(this1.i) > 0) {
        ohyea.i[1:length(this1.i),j] <- this1.i
      }
    }
    # change NA's back:
    ohyea[ohyea == "NA"] <- NA
    # remove rows of NA in ohyea.i:
    ohyea.i <- ohyea.i[rowSums(is.na(ohyea.i)) < ncol(ohyea.i),]
    
    if (return.conflix) {
      # rows (cols) of summary (index array) where there are conflicts (disagreeing names rather than target name+NA's)
      bloop <- which((rowSums(ohyea == "Bacteria", na.rm = TRUE) + rowSums(is.na(ohyea))) < (ncol(ohyea)-1))
      no.conflix.i <- unlist(ohyea.i[,-bloop]); no.conflix.i <- no.conflix.i[!is.na(no.conflix.i)] # indices to remove as bacteria...
      conflix.i <- unlist(ohyea.i[,bloop]); conflix.i <- conflix.i[!is.na(conflix.i)]
      if (length(conflix.i) > length(no.conflix.i)) {
        no.conflix.i <- c(no.conflix.i, rep(NA, times = length(conflix.i) - length(no.conflix.i)))
      } else if (length(conflix.i) < length(no.conflix.i)) {
        conflix.i <- c(conflix.i, rep(NA, times = length(no.conflix.i) - length(conflix.i)))
      }
      
      oh <- data.frame(n.conflict = c(length(conflix.i) - sum(is.na(conflix.i))), n.no.conflict = c(length(no.conflix.i) - sum(is.na(no.conflix.i))))
      yea <- data.frame(i.conflict = conflix.i, i.no.conflict = no.conflix.i)
      allout.list[[i]] <- list(oh, yea)
    } else {
      # store these arrays in a list within your big list
      allout.list[[i]] <- list(ohyea, ohyea.i)
    }
  }
  names(allout.list) <- taxnames
  
  return(allout.list)
}
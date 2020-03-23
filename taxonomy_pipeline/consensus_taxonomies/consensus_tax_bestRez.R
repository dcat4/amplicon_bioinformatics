consensus_tax_bestRez <- function(..., tablenames = c("bayes", "idtax"), ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"),
                          tiebreakz = "LCAlike"){
  nasum <- function(taxdf){
    x <- is.na(taxdf)
    ii <- rowSums(x)
    return(ii)
  }
  
  x <- list(...)
  notuz <- nrow(x[[1]])
  narray <- lapply(x, nasum) # the number of NA's in each row of each tax table
  narray <- matrix(unlist(narray), nrow=notuz, byrow=FALSE) # number na's in each taxonomy array
  bestrez <- apply(narray, 1, function(x) which(x == min(x))) # indices of the table (in x) that has the best resolution (least NA's) on a ASVxASV basis
  eh <- unlist(lapply(bestrez, length)) # looking for tie's where multiple taxonomies have = rez
  iTieBreak <- which(eh > 1)   # row indices where you'll need a tiebreaker (= resolution across tables)...
  # create a dataframe to store your final merged taxonomy:
  mergedTax <- data.frame(matrix(data=NA, nrow = notuz, ncol = length(ranknamez)))
  for(i in 1:length(tablenames)){
    bob <- which(unlist(lapply(bestrez, function(x) setequal(x,i)))) # indices where you can directly transfer taxonomies from table i in O.G. list
    fillerz <- x[[i]][bob,]
    mergedTax[bob,] <- fillerz
  }
  for (i in 1:length(iTieBreak)) {
    matchz <- lapply(x, function(x) x[iTieBreak[i],]) # all tax entries of a given tiebreaker index
    eh <- unique(matrix(unlist(matchz), nrow=length(matchz), byrow=TRUE)) # unique rows of tax assignments for this sequence
    # if there's only 1 unique batch of assignments for this row, all tax arrays agree.
    # in that case, add this unique assignment where it belongs in mergedTax:
    if (nrow(eh) == 1) {
      mergedTax[iTieBreak[i],] <- eh
      iTieBreak[i] <- NA # set to NA and remove later b/c no longer requires tiebreaking...
    } 
  }
  iTieBreak <- iTieBreak[!is.na(iTieBreak)]
  # this settles the tiebreakers based on input tiebreakz:
  if (length(grep("none", tiebreakz)) > 0) {
    # do nothing, no tiebreakerz
  } else if (identical("LCAlike", tiebreakz)) {
    # a loop that operates similar to LCA:
    for (i in 1:length(iTieBreak)) {
      tsi <- unlist(bestrez[iTieBreak[i]]) # indices of the tied taxonomy arrays
      xsub <- x[tsi]
      eh <- lapply(xsub, function(z) z[iTieBreak[i],]) # this ASV's taxonomies from each tied tax table
      eh <- lapply(eh, function(z) z[,!is.na(z)]) # remove NA ranks to get an accurate rank to assign to
      durp <-  xsub[[1]][iTieBreak[i],] # one taxonomy arrays assignments for comparison
      eh2 <- unlist(lapply(lapply(eh, intersect, durp), length)) # shortest entry will have probable LCA assignment...
      # there is an LCA, so assign it:
      lcai <- which(eh2 == min(eh2)) # index of tax array to use in LCA assignment
      if (length(lcai) > 1) {
        lcai <- min(lcai)
      }
      lcapathi <- min(eh2) # index of rank to assign to in LCA assignment
      if (lcapathi > 0) {
        mergedTax[iTieBreak[i],1:lcapathi] <- xsub[[lcai]][iTieBreak[i],1:lcapathi]
        iTieBreak[i] <- NA # set to NA and remove later b/c no longer requires tiebreaking...
      } else {
        # do nothing, LCA is at the root
        iTieBreak[i] <- NA # do remove this tiebreaker since you have technically assigned a root LCA
      }
    }
  } else {
    # this utilizes list input of tiebreakz described above
    for (j in 1:length(tiebreakz)) {
      pikl <- tiebreakz[[j]][1] # should be one of the tablenames or "LCAlike"
      # if it's LCAlike, do that routine:
      if (length(grep("LCAlike", pikl)) > 0) {
        # a loop that operates similar to LCA:
        for (i in 1:length(iTieBreak)) {
          tsi <- unlist(bestrez[iTieBreak[i]]) # indices of the tied taxonomy arrays
          xsub <- x[tsi]
          eh <- lapply(xsub, function(z) z[iTieBreak[i],]) # this ASV's taxonomies from each tied tax table
          eh <- lapply(eh, function(z) z[,!is.na(z)]) # remove NA ranks to get an accurate rank to assign to
          
          durp <-  xsub[[1]][iTieBreak[i],] # one taxonomy arrays assignments for comparison
          eh2 <- unlist(lapply(lapply(eh, intersect, durp), length)) # shortest entry will have probable LCA assignment...
          # there is an LCA, so assign it:
          lcai <- which(eh2 == min(eh2)) # index of tax array to use in LCA assignment
          if (length(lcai) > 1) {
            lcai <- min(lcai)
          }
          lcapathi <- min(eh2) # index of rank to assign to in LCA assignment
          
          if (lcapathi > 0) {
            mergedTax[iTieBreak[i],1:lcapathi] <- xsub[[lcai]][iTieBreak[i],1:lcapathi]
            iTieBreak[i] <- NA # set to NA and remove later b/c no longer requires tiebreaking...
          } else {
            # do nothing, LCA is at the root
            iTieBreak[i] <- NA # do remove this tiebreaker since you have technically assigned a root LCA
          }
        }
      } else {
        # if it's not LCAlike and not none, use the user-specified tax table:
        for (i in 1:length(iTieBreak)) {
          tsi <- unlist(bestrez[iTieBreak[i]]) # indices of the tied taxonomy arrays
          tnsub <- tablenames[tsi] # table names of tied tax arrays
          xsub <- x[tsi]
          eh <- which(tnsub %in% pikl) # index of tax array within xsub to use for this round of tiebreaking
          pn <- tiebreakz[[j]][2] # the taxonomic name to prioritize (or NA)
          # if eh is empty, save this i for another tiebreaking rule
          # if pn has no name just use the tiebreaking taxonomy
          # if pn has a name in it check for it in tiebreaking tax, use it if it's there, don't if it's not.
          if (length(eh) == 0) {
            # do nothing, you're tie-breaker doesn't apply here...
          } else if (length(eh) > 0 & is.na(pn)) {
            # assign the taxonomy with the tiebreaking array regardless of name:
            mergedTax[iTieBreak[i],] <- xsub[[eh]][iTieBreak[i],]
            iTieBreak[i] <- NA # set to NA and remove later b/c no longer requires tiebreaking...
          } else if (length(eh) > 0 & !is.na(pn)) {
            # check that this row of mergedTax contains the name specified..
            if(pn %in% xsub[[eh]][iTieBreak[i],]) {
              # if so do the assignment:
              mergedTax[iTieBreak[i],] <- xsub[[eh]][iTieBreak[i],]
              iTieBreak[i] <- NA # set to NA and remove later b/c no longer requires tiebreaking...
            } else {
              # if not do nothing
            }
          }
        }
      }
      iTieBreak <- iTieBreak[!is.na(iTieBreak)] # remove NAs from iTieBreak as those have had higher-priority tie-breakers applied...
    }
  }
  colnames(mergedTax) <- ranknamez
  y <- list(mergedTax, x, iTieBreak)
  return(y)
}


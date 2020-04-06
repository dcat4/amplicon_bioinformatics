consensus_tax_LCAlike <- function(..., tablenames = c("bayes", "idtax"), 
                                  ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species")) {
  x <- list(...)
  notuz <- nrow(x[[1]])
  mergedTax <- data.frame(matrix(data=NA, nrow = notuz, ncol = length(ranknamez)))
  for (i in 1:notuz) {
    eh <- lapply(x, function(z) z[i,]) # this ASV's taxonomies from each tied tax table
    eh <- lapply(eh, function(z) z[,!is.na(z)]) # remove NA ranks to get an accurate rank to assign to
    durp <-  x[[1]][i,] # one taxonomy arrays assignments for comparison
    eh2 <- unlist(lapply(lapply(eh, intersect, durp), length)) # shortest entry will have probable LCA assignment...
    lcai <- which(eh2 == min(eh2)) # index of tax array to use in LCA assignment
    if (length(lcai) > 1) {
      lcai <- min(lcai)
    }
    lcapathi <- min(eh2) # index of rank to assign to in LCA assignment
    if (lcapathi > 0) {
      mergedTax[i,1:lcapathi] <- x[[lcai]][i,1:lcapathi]
    } else {
      # do nothing, LCA is at the root
    }
  }
  colnames(mergedTax) <- ranknamez
  return(list(mergedTax,x))
}


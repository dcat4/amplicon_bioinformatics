# takes idtaxa outputs from silva (see dada2 tutorial) and formats them to match dada2's bayesian classifier outputs:

idtax2df_silva <- function(tt, boot = 0, 
                          ranksub = c("domain", "phylum", "class", "order", "family", "genus"), 
                          rubric = NULL, return.conf = FALSE){
  # adapted from dada2 tutorial
  tax <- t(sapply(tt, function(x) {
    m <- match(ranksub, x$rank)
    taxa <- x$taxon[m]
    taxa
  }))
  
  # adapted from dada2 tutorial
  conf <- t(sapply(tt, function(x) {
    m <- match(ranksub, x$rank)
    conf <- x$confidence[m]
    conf
  }))
  # this loop fills in NAs in between names (introduced by rank subsampling):
  # does so by mirroring dada2's filler names in bayesian-silva taxonomy:
  for (i in 1:nrow(tax)) {
    nana <- which(is.na(tax[i,]))
    bubu <- which(!is.na(tax[i,]))
    if (any(nana > bubu)) {
      for (j in 1:length(nana)) {
        ah <- nana[j] - bubu
        if (any(ah < 0)) {
          # you need to replace it with unclassified_ closest upstream name:
          eh <- which(ah == min(ah[ah > 0])) # this is the index in bubu of the closest name...
          suff <- substr(ranksub[nana[j]], start = 1, stop = 2) # suffix based on rank you're filling in
          tax[i, nana[j]] <- paste0(tax[i, bubu[eh]], "_", suff) 
          conf[i, nana[j]] <- conf[i, bubu[eh]] # propagate the confidence to these ranks too since you're not adding any real info...
        }
      }
    }
  }
  colnames(tax) <- ranksub
  colnames(conf) <- ranksub
  
  tax <- data.frame(tax, stringsAsFactors = FALSE)
  conf <- data.frame(conf, stringsAsFactors = FALSE)
  tax[conf < boot] <- NA # NA out below the boot threshold supplied
  rubdf <- data.frame(svN = names(rubric), ASV = as.character(rubric, use.names = FALSE), stringsAsFactors = FALSE)
  tax <- cbind(rubdf,tax)
  if (return.conf) {
    return(list(tax, conf))
  } else if (!return.conf) {
    return(tax)
  }
}
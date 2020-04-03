# this is a function that converts the output of the idtaxa algorithm as implemented in dada2 (see here: https://benjjneb.github.io/dada2/tutorial.html)
# outputs are a "taxon" object, which as far as I can tell is a list of list.
# regardless, the below will convert the taxon object returned by idtaxa (tt) into a dataframe, and will
# NA out any assignments with a confidence value less than boot (use 0 or no input [default=0] for no NA'ing)

# code written by Connie Liang (most) and Dylan Catlett (just polishing)

idtax2df <- function(tt, boot = 0, return.conf = "no"){
  taxonomy<-c()
  conf <- c()
  notu <- length(tt)
  for(j in 1:notu){
    taxonomy<-append(taxonomy,tt[[j]]$taxon)
    conf <- append(conf,tt[[j]]$confidence)
  }
  yydf <- data.frame(matrix(unlist(taxonomy), nrow=notu, byrow=TRUE), stringsAsFactors = FALSE)
  yydf <- yydf[,-1]
  confdf <- data.frame(matrix(unlist(conf), nrow=notu, byrow=TRUE))
  confdf <- confdf[,-1]
  yydf[confdf < boot] <- NA
  
  if (return.conf == "no") {
    return(yydf)
  } else if (return.conf == "yes") {
    return(list(yydf, confdf))
  }
}
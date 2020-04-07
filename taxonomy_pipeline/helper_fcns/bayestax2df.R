# converts the output of DADA2's assignTaxonomy to a dataframe you can plug into your pipeline
bayestax2df <- function(tt, boot = 0, rubric = NULL, return.conf = FALSE){
  taxonomy <- tt[[1]]
  conf <- tt[[2]]
  notu <- nrow(taxonomy)
  taxdf <- data.frame(taxonomy, stringsAsFactors = FALSE)
  
  confdf <- data.frame(conf, stringsAsFactors = FALSE)
  taxdf[confdf < boot] <- NA
  
  taxdf$ASV <- rownames(taxdf)
  rubdf <- data.frame(svN = names(rubric), ASV = as.character(rubric, use.names = FALSE), stringsAsFactors = FALSE)
  
  taxdf <- merge(taxdf, rubdf, by.x = "ASV", by.y = "ASV")
  taxdf <- taxdf[,c(ncol(taxdf), 1:(ncol(taxdf)-1))]
  if (return.conf) {
    return(list(taxdf, confdf))
  } else if (!return.conf) {
    return(taxdf)
  }
}
# converts the output of DADA2's assignTaxonomy to a dataframe you can plug into your pipeline
bayestax2df <- function(tt, boot = 0, rubric = NULL, return.conf = FALSE){
  taxonomy <- tt[[1]]
  conf <- tt[[2]]
  notu <- nrow(taxonomy)
  taxdf <- data.frame(taxonomy, stringsAsFactors = FALSE)
  
  confdf <- data.frame(conf)
  yydf[confdf < boot] <- NA
  
  yydf$ASV <- rownames(yydf)
  rubdf <- data.frame(svN = names(rubric), ASV = as.character(rubric, use.names = FALSE), stringsAsFactors = FALSE)
  
  yydf <- merge(yydf, rubdf, by.x = "ASV", by.y = "ASV")
  if (return.conf) {
    return(list(yydf, confdf))
  } else if (!return.conf) {
    return(yydf)
  }
}
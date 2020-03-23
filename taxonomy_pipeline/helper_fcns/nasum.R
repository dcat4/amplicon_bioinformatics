# this function returns a numeric vector of the number of NA values of each row of a dataframe

nasum <- function(taxdf){
  x <- is.na(taxdf)
  ii <- rowSums(x)
  return(ii)
}
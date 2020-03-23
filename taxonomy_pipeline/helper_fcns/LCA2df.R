# function to take LCA's output .csv and turn it into an R dataframe

LCA2df <- function(lcaer, rubric) {
  library("DECIPHER")
  cc <- apply(lcaer, MARGIN = 1, FUN = function(x) length(unlist(gregexpr(";", x))))
  nc <- max(cc)
  lcadf <- data.frame(matrix(NA, nrow = nrow(lcaer), ncol = nc))
  for (i in 1:nrow(lcaer)) {
    this1 <- lcaer[i,"V2"]
    nh <- unlist(gregexpr("No hits", this1))
    na <- unlist(gregexpr("Not assigned", this1))
    if (nh == -1 && na == -1) {
      # there's at least some hits...
      co <- unlist(gregexpr("cellular organisms", this1))
      ii <- unlist(gregexpr(";", this1))
      ii <- ii[ii > co]
      for (j in 1:(length(ii) - 1)) {
        lcadf[i,j] <- substr(this1,ii[j]+1,ii[j+1]-1)
      }
    } else {
      # leave it as all NA
    }
  }
  lcadf$svN <- la$V1 # append the SV numbers to your output dataframe
  rubdf <- data.frame(svN = names(rubric), ASV = as.character(rubric, use.names = FALSE))
  lcadf <- merge(lcadf, rubdf, by.x = "svN", by.y = "svN") # merging to add ASV sequences
  lcadf <- lcadf[,c(1,ncol(lcadf), 2:(ncol(lcadf)-1))] # rearranging columns -- svN, ASV, taxonomy....
  return(lcadf)
}
  

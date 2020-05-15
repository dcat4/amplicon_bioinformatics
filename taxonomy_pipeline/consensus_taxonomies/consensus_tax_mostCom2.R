consensus_tax_mostCom2 <- function(..., tablenames = c("bayes", "idtax"), ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"),
                                  tiebreakz = "LCAlike", count.na=FALSE, weights=rep(1, length(list(...)))) {
  x <- list(...)
  n.rows <- nrow(x[[1]])
  n.cols <- ncol(x[[1]])
  threshold <- 0.5
  n.dfs <- length(x)
  
  consensus.tax <- data.frame(matrix(ncol=n.cols, nrow=0, dimnames=list(NULL, names(x[[1]]))))
  
  # just iterate through each row and find the consensus of each row
  for (row in 1:n.rows) {
    # initialize the consensus row 
    c.row <- data.frame(matrix(rep(NA, n.cols), ncol=n.cols, nrow = 1, dimnames=list(NULL, names(consensus.tax))))
    c.row[, "svN"] <- x[[1]][row, "svN"]
    c.row[, "ASV"] <- x[[1]][row, "ASV"]
    for (col in rev(ranknamez)) {
      # collects the taxs from each data frame
      taxs <- vector()
      # corresponding vector to know which df the tax came from
      df.idx <- vector()
      for (i in 1:n.dfs) {
        df <- x[[i]]
        # weights are represented by the amount of repeated taxs are
        taxs <- c(taxs, rep(df[row, col], weights[i]))
        df.idx <- c(df.idx, rep(tablenames[i]), weights[i])
      }
      # create a frequency table to determine which one is the majority
      freq.df <- as.data.frame(table(taxs), stringsAsFactors=FALSE)
      # automatically, the frequency table excludes NA
      if (count.na) {
        # table parameter to count the NA's well
        freq.df <- as.data.frame(table(taxs, exclude=NULL), strinsgAsFactors=FALSE)
      }
      # if entires exist in the frequency table, determine the majority
      # else that means the NA's wasn't counted but will be set to NA as default
      if (nrow(freq.df) > 0) {
        # determine the proportion to compare to threshold by majority
        freq.df$prop <- freq.df$Freq / sum(freq.df$Freq)
        # get the one with the greatest proportion aka the majority
        # see if the proportion of the majority is at least the threshold
        max.prop <- max(freq.df$prop)
        if (max.prop >= threshold) {
          c.tax <- freq.df[which(freq.df$prop == max.prop), "taxs"]
          # if there are multiple taxs as the majority, we need to tie break it
          if (length(c.tax) > 1) {
            # see which data frame it is coming from
            # check if the data frame chosen is an option
            for (tax in c.tax) {
              idx <- which(tax %in% taxs)
              candidates <- df.idx[idx]
              if (iselement(tiebreaks, candidates)) {
                c.row[, col] <- tax 
                break
              }
            }
            # at this point just set it as NA
            c.row[, col] <- NA
          } else {
            c.row[, col] <- c.tax[1]
          }
        }
      }
    }
    # after iterating through the columns add the consensus row to the data frame
    consensus.tax <- rbind(consensus.tax, c.row)
  }
  
  return(consensus.tax)
}
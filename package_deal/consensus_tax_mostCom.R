# need to check on trueMajority argument -- doesn't look like it actually got written in...
# also frankenstein assignments are still possible when count.na = TRUE and not checked for in the error output

consensus_tax_mostCom <- function(..., tablenames = c("bayes", "idtax"), ranknamez = c("kingdom", "supergroup", "division","class","order","family","genus","species"),
                                   tiebreakz = "none", count.na=FALSE, trueMajority=FALSE, weights=rep(1, length(list(...)))) {
  library("dplyr")
  x <- list(...) # grab everything in a list structure 
  n.rows <- nrow(x[[1]]) # get the number of ASV's aka number of rows 
  n.cols <- ncol(x[[1]]) # get the number of columns to recreate consensus taxonomy table 
  threshold <- 0.5 # threshold preset
  n.dfs <- length(x) # getting the number of taxonomy tables inputted 
  
  # set up empty consensus taxonomy table to add by each row
  consensus.tax <- data.frame(matrix(ncol=n.cols, nrow=0, dimnames=list(NULL, names(x[[1]]))))
  
  # determine if tiebreaks were specified
  tiebreaker <- NA
  if (tiebreakz != "none") {
    # convert list of tiebreakers to a dataframe t
    tiebreaker <- data.frame(matrix(unlist(tiebreakz), nrow=length(tiebreakz), byrow=T, dimnames=list(NULL, c("table", "tax"))),stringsAsFactors=FALSE)
    # set priority numbers
    tiebreaker$priority <- as.numeric(rownames(tiebreaker))
  }
  
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
        df.idx <- c(df.idx, rep(tablenames[i], weights[i]))
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
        if (!trueMajority || max.prop >= threshold) {
          c.tax <- as.character(freq.df[which(freq.df$prop == max.prop), "taxs"])
          # if there are multiple taxs as the majority, we need to tie break it
          if (length(c.tax) > 1) {
            # see which data frame it is coming from
            # check if the data frame chosen is an option
            
            if (!is.na(tiebreaker)) {
              # create a data frame with taxs and tablename it came from
              pairs <- vector(mode="list", length=length(c.tax))
              curr.idx <- 1
              for (i in 1:length(c.tax)) {
                # find the corresponding tablename of tax
                tax <- c.tax[i]
                idx <- match(tax, taxs)
                if (is.na(tax)) {
                  tax <- "na"
                }
                table <- df.idx[idx]
                # append result to list
                pairs[[curr.idx]] <- c(table, tax)
                curr.idx <- curr.idx + 1
              }
              # data frame of ties
              ties <- data.frame(matrix(unlist(pairs), nrow=length(pairs), byrow=T, dimnames=list(NULL, c("table", "tax"))),stringsAsFactors=FALSE)
              
              # first assign exact matches
              exact <- merge(ties, tiebreaker, by=c("table","tax"), all.x=TRUE)
              
              # go back and consider taxnames only with NA
              all <- merge(exact, tiebreaker[which(is.na(tiebreaker[,"tax"])), c("table","priority")], by="table", all.x=TRUE)
              
              # resolve priorities
              all$priority <- coalesce(all$priority.x, all$priority.y)
              all$priority.x <- NULL
              all$priority.y <- NULL
              
              # sort by priority to get tiebreaker at the top of the data frame
              sorted <- all[order(all$priority), ]
              
              # assign top row as consensus
              if (is.na(sorted[1, "priority"])) {
                c.row[, col] <- NA
              } else {
                c.row[, col] <- sorted[1, "tax"]
              }
            } else {
              # at this point just set it as NA
              c.row[, col] <- NA
            }
          } else {
            c.row[, col] <- c.tax[1]
          }
        }
      }
    }
    # after iterating through the columns add the consensus row to the data frame
    consensus.tax <- rbind(consensus.tax, c.row)
  }
  
  # convert all NA strings to actual NA's
  make.true.NA <- function(x) if(is.character(x)||is.factor(x)){
    is.na(x) <- x=="NA"; x} else {
      x}
  
  df[] <- lapply(consensus.tax, make.true.NA)
  df <- data.frame(lapply(df, as.character), stringsAsFactors=FALSE)
  # dylan's addition - checking and warning for non-optimal assignments:
  # below returns TRUE where an NA is found in the middle of a heirarchical assignment for a particular ASV
  qcer <- function(y) {
    eh <- is.na(df)
    ina <- apply(eh, MARGIN = 1, FUN = function(x) min(which(x)))
    nna <- apply(eh, MARGIN = 1, FUN = function(x) max(which(!x)))
    return(any(ina < nna))
  }
  
  qc1 <- qcer(df)
  if (qc1) {
    stop(c("Non-optimal taxonomic assignments detected in consensus. \nI advise re-computing with 'count.na = TRUE'"))
  }
  return()
}
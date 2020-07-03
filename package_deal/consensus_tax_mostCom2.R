# this is a re-write of Kevin's mostCom after finding the frankenstein assignments are possible for basically all of parameter spaces available
# I think it should work now. Need to add in the Frankenstein checker, and then do some more tests.

# The only parameters that (I think) will not work are:
# 1. count.na = FALSE -- I think it will always count NA's ?
# 2. tie-breaking where you want to specify the a tax name with table ; i think specifying a table will still work.
# 3. test the count.na = FALSE

consensus_tax_mostCom2 <- function(..., tablenames = c(), ranknamez = c("kingdom", "supergroup", "division","class","order","family","genus","species"),
                                   tiebreakz = "none", count.na=TRUE, trueMajority=TRUE, threshold = 0.5, weights=c()) {
  library("dplyr")
  x <- list(...) # grab everything in a list structure 
  n.rows <- nrow(x[[1]]) # get the number of ASV's aka number of rows 
  n.cols <- ncol(x[[1]]) # get the number of columns to recreate consensus taxonomy table 
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
    # message("row = ", toString(row))
    c.row <- data.frame(matrix(rep(NA, n.cols), ncol=n.cols, nrow = 1, dimnames=list(NULL, names(consensus.tax))))
    c.row[, "svN"] <- x[[1]][row, "svN"]
    c.row[, "ASV"] <- x[[1]][row, "ASV"]
    # create an alltax dataframe to track heirarchical assignments:
    alltax <- data.frame(matrix(NA, nrow = n.dfs, ncol = n.cols, dimnames=list(names(x), names(consensus.tax))), stringsAsFactors = FALSE)
    for (i in 1:n.dfs) {
      df <- x[[i]]
      alltax[i , ] <- df[row , ]
    }
    for (col in ranknamez) {
      
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
        freq.df <- as.data.frame(table(taxs, exclude=NULL), stringsAsFactors=FALSE)
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
              pairs <- data.frame()
              for (i in 1:length(c.tax)) {
                # find the corresponding tablename of tax
                tax <- c.tax[i]
                idx <- taxs %in% tax 
                tbl <- df.idx[idx]
                pairs <- rbind(pairs, cbind(tbl, rep(tax, times = length(tbl))))
              }
              # data frame of ties
              colnames(pairs) <- c("table","tax")
              # first assign exact matches
              exact <- merge(pairs, tiebreaker, by=c("table","tax"), all.x=TRUE)
              
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
            c.row[, col] <- c.tax
          }
        }
      }
      # if the assignment was NA, break (all downstream ranks are already NA)
      if (is.na(c.row[, col])) {
        break
      }
    }
    c.row <- data.frame(lapply(c.row, as.character), stringsAsFactors=FALSE)
    
    # The below ensures that the consensus is derived from the inputs and not 'frankensteined'
    tmp1 <- c.row
    tmp2 <- alltax
    checker <- dplyr::intersect(tmp1, tmp2)
    while (nrow(checker) == 0){ 
      tmp1 <- tmp1[, -ncol(tmp1)]
      tmp2 <- tmp2[, -ncol(tmp2)]
      checker <- dplyr::intersect(tmp1, tmp2)
    }
    if (ncol(checker) < n.cols) {
      c.row[, (ncol(checker)+1):n.cols] <- NA
    }
    
    # after iterating through the columns add the consensus row to the data frame
    consensus.tax <- rbind(consensus.tax, c.row)
  }
  
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
    message("Frankenstein assignments detected. You should modify your ensemble parameters and try again.")
  }
  tmp <- c(x, list(consensus.tax))
  names(tmp) <- c(tablenames, "ensemble")
  check4frankenstein(tmp, ranknames = ranknamez)
  return(df)
}

check4frankenstein <- function(tbl.list, ranknames = c("kingdom","supergroup","division","class","order","family","genus","species")) {
  # pull out ensemble and the output of each individual algorithm:
  ee <- tbl.list[["ensemble"]]
  tt <- tbl.list[names(tbl.list) != "ensemble"]
  tt <- bind_rows(tt)
  # extract unique taxonomic paths from the ensemble and all individual tables:
  uee <- unique(ee, MARGIN = 1)
  tl <- unique(tt, MARGIN = 1)
  
  uee <- ee
  tl <- tt
  nchex <- nrow(uee)
  for (row in 1:nchex) {
    checker <- dplyr::intersect(uee[row , ], tl)
    if (nrow(checker) == 0){
      ee <- uee[row , ]
      tmp <- tl
      tmp <- tmp[, 1:max(which(!is.na(ee)))]
      ee <- ee[, 1:max(which(!is.na(ee)))]
      check2 <- dplyr::intersect(ee, tmp)
      if (nrow(check2) == 0){
        # this is a frankensteined assignment. throw an error and break
        message("Frankenstein assignments detected. You should modify your ensemble parameters and try again.")
      }
    }
  }
  # if you make it here, there are no frankenstein assignments. print that and you're done
  message("no Frankenstein assignments detected")
}

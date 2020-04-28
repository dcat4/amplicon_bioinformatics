consensus_tax_mostCom <- function(..., tablenames = c("bayes", "idtax"), ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"),
                                  tiebreakz = "LCAlike", count.na=FALSE, weights=rep(1, length(list(...)))) {
  x <- list(...)
  consensus.tax <- data.frame(matrix(ncol=(ncol(x[[1]])), nrow=0, dimnames=list(NULL, names(x[[1]]))))
  n.rows <- nrow(x[[1]])
  n.cols <- ncol(x[[1]])
  threshold <- 0.5
  n.dfs <- length(x)
  ties <- vector()
  
  # get the rows that are all NA's and pop them into the consensus.tax since they're all NA's
  # think of ways to validate the mapping of consensus
  
  for (row in 1:n.rows) {
    placed <- FALSE
    for (col in rev(3:n.cols)) {
      taxs <- vector()
      for (i in 1:n.dfs) {
        df <- x[[i]]
        taxs <- c(taxs, rep(df[row, col], weights[i]))
      }
      freq.df <- as.data.frame(table(taxs), stringsAsFactors=FALSE)
      # add NA freq here by counting NA in taxs
      if (count.na) {
        freq.df <- as.data.frame(table(taxs, exclude=NULL), stringsAsFactors=FALSE)
        # if the table only contains NA move to the next column
        if (nrow(freq.df) == 1 & (is.na(freq.df[,1]))) {
          next
        }
      }
      if (nrow(freq.df) > 0) {
        # multiply the weight to the frequency and divide it by the normalized 
        
        # example of weighing (default vector of 1's)
        # table1 -> Bacteria      table2 -> Bacteria      table3 -> Eukaryota 
        # weighing scheme 1 1 2
        
        # freq -> 2 Euk and 2 Bacteria 
        
        
        freq.df$prop <- freq.df$Freq / sum(freq.df$Freq)
        max.prop <- max(freq.df$prop)
        if (max.prop >= threshold) {
          c.tax <- freq.df[which(freq.df$prop == max.prop), "taxs"]
          if (length(c.tax) > 1) {
            ties <- c(ties, row)
            placed <- TRUE
            break
          }
          else {
            # add NA downstream of the column found 
            df <- x[[match(c.tax, taxs)]]
            
            matched.row <- data.frame(matrix(rep(NA, n.cols), ncol = n.cols, nrow = 1))
            colnames(matched.row) <- colnames(df)
            matched.row[1:col] <- df[row, ][1:col]
            
            consensus.tax <- rbind(consensus.tax, matched.row)
            placed <- TRUE
            break
          }
        }
      }
    }
    if(!placed) {
      ties <- c(ties, row)
    }
  }
  
  # tie breaking
  # tie for the highest proportion
  # take the tables -> subset it and apply the tie breaking rules to the subset
  
  for (i in 1:length(ties)) {
    matchz <- lapply(x, function(x) x[ties[i], ])
    eh <- unique(matrix(unlist(matchz), nrow=length(matchz), byrow=TRUE))
    if (nrow(eh) == 1) {
      consensus.tax <- rbind(consensus.tax, eh)
      ties[i] <- NA
    }
  }
  ties <- ties[!is.na(ties)]
  if (length(grep("none", tiebreakz)) > 0) {

  } else if (identical("LCAlike", tiebreakz)) {
    for (i in 1:length(ties)) {
      eh <- lapply(x, function(z) z[ties[i],])
      eh <- lapply(eh, function(z) z[, !is.na(z)])
      durp <- x[[1]][ties[i],]
      eh2 <- unlist(lapply(lapply(eh, intersect, durp), length))
      lcai <- which(eh2 == min(eh2))
      if (length(lcai) > 1) {
        lcai <- min(lcai)
      }
      if (lcapathi > 0) {
        matched.row <- data.frame(matrix(rep(NA, n.cols), ncol = n.cols, nrow = 1))
        colnames(matched.row) <- colnames(df)
        matched.row[,1:lcapathi] <- x[[lcai]][ties[i], 1:lcapathi]
        consensus.tax <- rbind(consensus.tax, matched.row)
        ties[i] <- NA
      } else {
        ties[i] <- NA
      }
    }
  } else {
    for (j in 1:length(tiebreakz)) {
      pikl <- tiebreakz[[j]][1]
      if (length(grep("LCAlike", pikl)) > 0) {
        for (i in 1:length(ties)) {
          eh <- lapply(x, function(z) z[ties[i],])
          eh <- lapply(eh, function(z) z[,!is.na(z)])
          durp <- x[[1]][ties[i],]
          eh2 <- unlist(lapply(lapply(eh, intersect, durp), length))
          lcai <- which(eh2 == min(eh2))
          if (length(lcai) > 1) {
            lcai <- min(lcai)
          }
          lcapathi <- min(eh2)
          if (lcapathi > 0) {
            matched.row <- data.frame(matrix(rep(NA, n.cols), ncol = n.cols, nrow = 1))
            colnames(matched.row) <- colnames(df)
            matched.row[,1:lcapathi] <- x[[lcai]][ties[i], 1:lcapathi]
            consensus.tax <- rbind(consensus.tax, matched.row)
            ties[i] <- NA
          } else {
            ties[i] <- NA
          }
        }
      } else {
        # table names ?
      }
      ties <- ties[!is.na(ties)]
    }
  }
  
  # add the tie breaking stuff -> user specified rules
  # weighing option for each data frame for proportion calculation (numeric vector in order of df)
  # option to treat NA's as a name
  
  return(list(consensus.tax, ties))
  
  # NOTES
  # majority rule system
  # use the assignment that is in agreement in the most taxonomy tables
  # row and rank wise basis majority rule system if the majority of the input taxonomy
  
  # try to avoid for loops (lapply?)
  
  # assume homoegenous ranking convention -> each input follow the same ranking
  
  # have an arbitrary number of input of taxonomy tables where they are all formatted the same way
  # row 1 corresponds to row 1 of another table 
  
  # assuming that you have three input tables...
  # for each table in table
  # row 1 would be ASV1 and would have a series of assignments 
  # for each of those name
  # if this name is shared by two or more of the tables, retain the name (majority)
  
  # for each row, start at the highest column first and iterate backwards to first column
  # check for equivalency -> if equal, pick one row and assign to the output
  # if two of them are equal -> assign either 
  # if > 50% of the matching  -> assign it to the upstream name 
  # < 50% -> use tie breaking scheme in consensus -> copy/paste 
  # store those row indicies and loop through a tie breaking loop 
  
  # consider A = tax, B = NA, C = NA 
  # have another input if NA is considered as a NA
  # -> if that input is true, value is NA
  # have some input where you have one option as NA as standard name -> assign upstream names
  
  # start treating NA's as names
  # do it by row
  # store it for the tie breakers if anything downstream doesn't match or use NA
}
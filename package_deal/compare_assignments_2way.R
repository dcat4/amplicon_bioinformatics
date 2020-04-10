# written by D Catlett
# code contributed by C Liang
# will be used as an automated function to compare taxonomic assignments made by 2 taxonomy tables.

compare_assignments_2way <- function(..., pltfilez = "none",
                                     tablenames = c("bayes", "idtax"), 
                                     ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species")) {
  library("ggplot2")
  nasum <- function(taxdf) {
    notuz <- nrow(taxdf)
    x <- is.na(taxdf)
    ii <- rowSums(x)
    return(ii)
  }
  x <- list(...)
  notuz <- nrow(x[[1]]) # number of ASVs/OTUs/rows in each tax table
  table1 <- x[[1]]
  table2 <- x[[2]]
  # you can find rows in which one taxonomy table has more NA's than the other by lapply-ing nasum:
  narray <- lapply(x, nasum)
  # vectors for storing indices of various comparison outcomes:
  perf <- c() # perfect agreement
  sameN.mo1R.i <- c() # same names, table 1 more resolved
  sameN.mo2R.i <- c() # same names, table 2 more resolved
  sameN.mo1R.r <- c() # same names, table 1 more resolved -- stores rank of less resolved table2
  sameN.mo2R.r <- c() # same names, table 2 more resolved -- stores rank of less resolved table1
  diffN <- c() # no names in common
  same2r.i <- c() # same names to a particular rank, after that disagreement -- index vector
  same2r.r <- c() # same names to a particular rank, after that there's disagreement -- first rank at which they disagree
  # bothna <- c() # for storing indices of rows where both tables are entirely unresolved
  for (i in 1:notuz) {
    # extract ith ASv from each table and remove NA's:
    t1 <- table1[i,]
    t2 <- table2[i,]
    t1 <- t1[,!is.na(t1)] 
    t2 <- t2[,!is.na(t2)]
    if (identical(t1,t2) || (length(setdiff(t1,t2)) == 0 && length(setdiff(t2,t1)) == 0)) {
      
      perf <- append(perf,i) # perfect agreement
      # THIS WILL INCLUDE ASSIGNMENTS THAT ARE ALL NA IN BOTH ARRAYS!!
    } else if (length(intersect(t1,t2)) == 0) {
      if (length(t1) == 0 && length(t2) == 0) {
        feck
        # this should never happen, so break if it does.
      } else if (length(t1) > 0 && length(t2) > 0) {
        # this means they have no names in common and neither were entirely unassigned
        # that's a perfect disagreement:
        diffN <- append(diffN, i)
      } else if (length(t1) > 0 && length(t2) == 0) {
        # this means t2 is entirely unassigned, but t1 has something:
        sameN.mo1R.i <- append(sameN.mo1R.i,i)
        sameN.mo1R.r <- append(sameN.mo1R.r,length(t2))
      } else if (length(t1) == 0 && length(t2) > 0) {
        # this means t1 is entirely unassigned, but t2 has something:
        sameN.mo2R.i <- append(sameN.mo2R.i,i)
        sameN.mo2R.r <- append(sameN.mo2R.r,length(t1))
      }

    } else if (length(intersect(t1,t2)) > 0) {
      # this means they have >= 1 name agreement, potentially varying resolution
      if (narray[[1]][i] > narray[[2]][i]) {
        # table 2 has better resolution...
        bloop <- length(setdiff(t1,t2)) # number of different, non-NA names 
        # if less-resolved table (t1 here) is first and all names it has are found in more resolved table, bloop will be 0.
        if (bloop == 0){
          # there are no disagreeing names, so populate accordingly:
          sameN.mo2R.i <- append(sameN.mo2R.i,i)
          sameN.mo2R.r <- append(sameN.mo2R.r,length(t1))
        } else if (bloop > 0) {
          # there are bloop disagreeing names - use that to inform the rank to which they agree...
          same2r.i <- append(same2r.i,i)
          same2r.r <- append(same2r.r, length(t1) - bloop)
        }
      } else if (narray[[1]][i] < narray[[2]][i]) {
        # table 1 has better resolution...
        bloop <- length(setdiff(t2,t1)) # number of different, non-NA names
        if (bloop == 0){
          # there are no disagreeing names, so populate accordingly:
          sameN.mo1R.i <- append(sameN.mo1R.i,i)
          sameN.mo1R.r <- append(sameN.mo1R.r,length(t2))
        } else if (bloop > 0) {
          # there are bloop disagreeing names - use that to inform the rank to which they agree...
          same2r.i <- append(same2r.i,i)
          same2r.r <- append(same2r.r, length(t2) - bloop)
        }
      } else if (narray[[1]][i] == narray[[2]][i]) {
        # they have = resolution, so they must have at least 1 different name
        same2r.r <- append(same2r.r, length(intersect(t1,t2)))
        same2r.i <- append(same2r.i,i)
      } 
    }
    # still in the for-loop here..
  }
  # wrap up your outputs...
  # concatenate index vectors into a single dataframe:
  nr <- max(c(length(same2r.i),length(sameN.mo1R.i),length(sameN.mo2R.i),length(perf),length(diffN)))
  indexDF <- data.frame(matrix(NA, nrow = nr, ncol = 8))
  colnames(indexDF) <- c("Agree", "Disagree", "Agree.T1.more.rez", "T1.more.rez.t2rank", "Agree.T2.more.rez", "T2.more.rez.t1rank", "pAgree.i", "pAgree.rank")
  if (length(perf) > 0){
    indexDF[1:length(perf),"Agree"] <- perf
  }
  if (length(diffN) > 0) {
    indexDF[1:length(diffN),"Disagree"] <- diffN
  }
  if (length(sameN.mo1R.i) > 0) {
    indexDF[1:length(sameN.mo1R.i),"Agree.T1.more.rez"] <- sameN.mo1R.i
    indexDF[1:length(sameN.mo1R.r),"T1.more.rez.t2rank"] <- sameN.mo1R.r
  }
  if (length(sameN.mo2R.i) > 0) {
    indexDF[1:length(sameN.mo2R.i),"Agree.T2.more.rez"] <- sameN.mo2R.i
    indexDF[1:length(sameN.mo2R.r),"T2.more.rez.t1rank"] <- sameN.mo2R.r
  }
  if (length(same2r.i) > 0) {
    indexDF[1:length(same2r.i),"pAgree.i"] <- same2r.i
  }
  if (length(same2r.r) > 0) {
    indexDF[1:length(same2r.r),"pAgree.rank"] <- same2r.r
  }
  plotDF <- data.frame(matrix(NA, nrow = 5, ncol = 2))
  colnames(plotDF) <- c("comp.val", "nOTU")
  plotDF[,"nOTU"] <- c(length(perf),length(diffN),length(sameN.mo1R.i),length(sameN.mo2R.i),length(same2r.i))
  plotDF[,"propOTU"] <- plotDF[,"nOTU"] / notuz
  plotDF[,"comp.val"] <- c("Total Agreement", "Total Disagreement", 
                           paste0("Names Agree,\n", tablenames[1], " - Higher Rez"), paste0("Names Agree,\n", tablenames[2], " - Higher Rez"), 
                           "Partial Agreement")
  # plotting with ggplot2:
  library("ggplot2")
  # proportional number of ASVs
  p1 <- ggplot(plotDF, aes(x = comp.val, y = propOTU)) + 
    geom_bar(stat="identity", color = "black", position=position_dodge()) + 
    labs(x = "", y = "Proportion of ASVs") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_text(size = 12, face="bold"),
          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12, face="bold"),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "white"),
          axis.line = element_line(size = 0.5, linetype = "solid", colour = "black")) +
    ggtitle(paste0(tablenames[1], " vs. ",  tablenames[2]))
  
  # absolute number of ASVs:
  p2 <- ggplot(plotDF, aes(x = comp.val, y = nOTU)) + 
    geom_bar(stat="identity", color = "black", position=position_dodge()) + 
    labs(x = "", y = "N ASVs") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_text(size = 12, face="bold"),
          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12, face="bold"),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "white"),
          axis.line = element_line(size = 0.5, linetype = "solid", colour = "black")) +
    ggtitle(paste0(tablenames[1], " vs. ",  tablenames[2]))
  
  if (pltfilez == "none"){
    # don't save anything
  } else {
    ggsave(filename = pltfile, plot = p1, device = "pdf")
    ggsave(filename = pltfile, plot = p2, device = "pdf")
  }
  return(list(indexDF, plotDF, p1, p2))
}

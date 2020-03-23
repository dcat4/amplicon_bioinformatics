# under construction...

# expanding on the 2-way comparisons w/ 3-way comparisons:
# for now litereally copy/pasted from 2way with minimal tweaks.

# strategizing:
# maybe wrap the 2way function on pairwise comparisons and that would trim the outcomes some...

# Categories:
# 1. 3x perfect agreement -- can do w/ 2way and intersection of agree
#   a. 2x perfect + 1 partial
#   b. 2x perfect + 1 rez
#   c. 2x perfect + 1 dis
# 2. 3x partial agreements -- can prob do w/ 2way and intersection of partial agreements...
#   a. should include 
# 3. resolution variability -- can prob do w/ 2way at least partly
# 4. perfect disagreement -- can't do by wrapping 2way I don't think

compare_assignments_3way <- function(..., pltfilez = c("prop_2wayplt.pdf", "abs_2wayplt.pdf"),
                                     tablenames = c("bayes", "idtax"), 
                                     ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species")) {
  library("ggplot2")
  
  x <- list(...)
  notuz <- nrow(x[[1]]) # number of ASVs/OTUs/rows in each tax table
  table1 <- x[[1]]
  table2 <- x[[2]]
  table3 <- x[[3]]
  
  # read in the 2-way comparison function:
  source("~/Documents/R/desktop_ampData_processing/connie_taxonomy_stuff_Mar2020/18sV9_amplicon_sequencing/taxonomy_pipeline/tax_table_comparisons/compare_assignments_2way.R")
  # we'll compare these 2 to start:
  tblnam <- c("t1", "t2")
  c2.12 <- compare_assignments_2way(table1, table2, tablenames = tblnam, pltfilez = "none")
  c2.13 <- compare_assignments_2way(table1, table3, tablenames = tblnam, pltfilez = "none")
  c2.23 <- compare_assignments_2way(table2, table3, tablenames = tblnam, pltfilez = "none")
  idf.12 <- c2.12[[1]]
  idf.13 <- c2.13[[1]]
  idf.23 <- c2.23[[1]]

  # get perfect indices by intersecting agreement indices in the 1-2 and 2-3 comparisons:
  perf3 <- intersect(idf.12$Agree, idf.23$Agree)
  
  return(list(x, idf.12, idf.13, idf.23))
  # LEFT OFF HERE IN BUILDING...
  
  # vectors for storing indices of various comparison outcomes:
  sameN.mo1R <- c() # same names, table 1 more resolved
  sameN.mo2R <- c() # same names, table 2 more resolved
  sameN.mo3R <- c() # same names, table 3 more resolved
  
  diffN <- c() # both resolved, no names in common
  
  same2r.i <- c() # same names to a particular rank, after that disagreement -- index vector
  same2r.r <- c() # same names to a particular rank, after that there's disagreement -- first rank at which they disagree
  
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
        sameN.mo1R <- append(sameN.mo1R,i)
      } else if (length(t1) == 0 && length(t2) > 0) {
        # this means t1 is entirely unassigned, but t2 has something:
        sameN.mo2R <- append(sameN.mo2R,i)
      }
      
    } else if (length(intersect(t1,t2)) > 0) {
      # this means they have >= 1 name agreement, potentially varying resolution
      if (narray[[1]][i] > narray[[2]][i]) {
        # table 2 has better resolution...
        bloop <- length(setdiff(t1,t2)) # number of different, non-NA names 
        # if less-resolved table (t1 here) is first and all names it has are found in more resolved table, bloop will be 0.
        if (bloop == 0){
          # there are no disagreeing names, so populate accordingly:
          sameN.mo2R <- append(sameN.mo2R,i)
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
          sameN.mo1R <- append(sameN.mo1R,i)
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
  nr <- max(c(length(same2r.i),length(sameN.mo1R),length(sameN.mo2R),length(perf),length(diffN)))
  indexDF <- data.frame(matrix(NA, nrow = nr, ncol = 6))
  colnames(indexDF) <- c("Agree", "Disagree", "Agree.T1.more.rez", "Agree.T2.more.rez", "pAgree.i", "pAgree.rank")
  if (length(perf) > 0){
    indexDF[1:length(perf),"Agree"] <- perf
  }
  if (length(diffN) > 0) {
    indexDF[1:length(diffN),"Disagree"] <- diffN
  }
  if (length(sameN.mo1R) > 0) {
    indexDF[1:length(sameN.mo1R),"Agree.T1.more.rez"] <- sameN.mo1R
  }
  if (length(sameN.mo2R) > 0) {
    indexDF[1:length(sameN.mo2R),"Agree.T2.more.rez"] <- sameN.mo2R
  }
  if (length(same2r.i) > 0) {
    indexDF[1:length(same2r.i),"pAgree.i"] <- same2r.i
  }
  if (length(same2r.r) > 0) {
    indexDF[1:length(same2r.r),"pAgree.rank"] <- same2r.r
  }
  plotDF <- data.frame(matrix(NA, nrow = 5, ncol = 2))
  colnames(plotDF) <- c("comp.val", "nOTU")
  plotDF[,"nOTU"] <- c(length(perf),length(diffN),length(sameN.mo1R),length(sameN.mo2R),length(same2r.i))
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
  
  if (length(grep(pltfilez, "none")) == 1){
    # don't save anything
  } else {
    ggsave(filename = pltfile, plot = p1, device = "pdf")
    ggsave(filename = pltfile, plot = p2, device = "pdf")
  }
  
  return(list(indexDF, plotDF, p1, p2))
}

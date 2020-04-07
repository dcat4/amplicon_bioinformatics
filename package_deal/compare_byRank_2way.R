# this function compares two taxonomy tables in a rank-wise fashion. 

compare_byRank_2way <- function(table1, table2,
                                pltfilez = "none",
                                tablenames = c("bayes", "idtax"), 
                                ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species")) {
  notuz <- nrow(table1) # number of ASVs/OTUs/rows in each tax table
  # initialize a list for storing index vectors:
  allofit <- list()
  # you have to apply the below to each column, potentially in a loop
  for (i in 1:length(ranknamez)) {
    i.bothna <- intersect(which(is.na(table1[,i])), which(is.na(table2[,i])))
    i.t1N.t2na <- intersect(which(!is.na(table1[,i])), which(is.na(table2[,i])))
    i.t1na.t2N <- intersect(which(is.na(table1[,i])), which(!is.na(table2[,i])))
    i.bothN <- intersect(which(!is.na(table1[,i])), which(!is.na(table2[,i])))
    i.bothN.a <- intersect(i.bothN, which(table1[,i] == table2[,i]))
    i.bothN.d <- intersect(i.bothN, which(table1[,i] != table2[,i]))
    # prob merge these into a dataframe and make each entry in index list correspond to a rank...
    nr <- max(c(length(i.bothna),length(i.t1N.t2na),length(i.t1na.t2N),length(i.bothN.a),length(i.bothN.d)))
    indexDF <-  data.frame(matrix(NA, nrow = nr, ncol = 5))
    colnames(indexDF) <- c("both.na", "t1named.t2NA", "t2named.t1NA", "same.name", "diff.name")
    if (length(i.bothna) > 0){
      indexDF[1:length(i.bothna), "both.na"] <- i.bothna
    }
    if (length(i.t1N.t2na) > 0){
      indexDF[1:length(i.t1N.t2na), "t1named.t2NA"] <- i.t1N.t2na
    }
    if (length(i.t1na.t2N) > 0){
      indexDF[1:length(i.t1na.t2N), "t2named.t1NA"] <- i.t1na.t2N
    }
    if (length(i.bothN.a) > 0){
      indexDF[1:length(i.bothN.a), "same.name"] <- i.bothN.a
    }
    if (length(i.bothN.d) > 0){
      indexDF[1:length(i.bothN.d), "diff.name"] <- i.bothN.d
    }
    pp <- apply(indexDF, MARGIN = 2, function(x) length(which(!is.na(x)))) # one rank's worth of numbers you can plot
    if (i == 1) {
      plotDF <- data.frame(comp = names(pp), count = pp, rank = rep(ranknamez[i], times = length(pp)))
    } else {
      plotDF <- rbind(plotDF, data.frame(comp = names(pp), count = pp, rank = rep(ranknamez[i], times = length(pp))))
    }
    allofit[[i]] <- indexDF
  }
  # outside the loop now
  
  # move to plotting:
  plotDF[,"count"] <- plotDF[,"count"] / notuz
  
  # plotting with ggplot2:
  library("ggplot2")
  # comp on x, rank on color
  p1 <- ggplot(plotDF, aes(x = rank, y = count, fill = comp)) +
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
    scale_fill_discrete(name = "Comparison", labels = c("Both NA", "Different Name", "Same Name", paste0(tablenames[1], " named,\n", tablenames[2], " NA"), paste0(tablenames[2], " named,\n", tablenames[1], " NA"))) +
    ggtitle(paste0(tablenames[1], " vs. ",  tablenames[2]))
  
  # comp on color, rank on x
  p2 <- ggplot(plotDF, aes(x = comp, y = count, fill = rank)) + 
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
    scale_fill_discrete(name = "Rank") +
    scale_x_discrete(labels = c(both.na = "Both NA", diff.name = "Different Name", same.name = "Same Name", t1named.t2NA = paste0(tablenames[1], " named,\n", tablenames[2], " NA"), t2named.t1NA = paste0(tablenames[2], " named,\n", tablenames[1], " NA")))
    ggtitle(paste0(tablenames[1], " vs. ",  tablenames[2]))
  
  if (length(grep(pltfilez, "none")) == 1){
    # don't save anything
  } else {
    ggsave(filename = pltfile, plot = p1, device = "pdf")
    ggsave(filename = pltfile, plot = p2, device = "pdf")
  }
  return(list(allofit, plotDF, p1, p2))
}
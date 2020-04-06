analyze_traitmap_byTrait <- function(map.result, trait.name, otu.table = "none", plotfilez = "none") {
  library("stringr")
  library("ggplot2")
  library("reshape2")
  tvec <- map.result[,trait.name]
  tvec.na <- is.na(tvec) # where the name was mapped directly to NA, and nothing else
  tvec.containa <- str_which(tvec, "NA") # where the hit was NA in some instances but something else in others...
  nhitz <- str_count(tvec, ";")
  # nhitz will be 0 if there was 1 unique value for the trait
  # otherwise it will == the number of unique values of this trait
  i.am <- sort(c(which(tvec.na), which(nhitz > 0))) # ambiguous
  i.unam <- which(nhitz == 0) # unambiguous
  
  bloop <- unique(tvec[-c(tvec.containa, which(tvec.na))]) # unique mapping outcomes excluding entries containing NA
  # ncol below should be length(bloop) + 1 (outcomes w/ NA) + 2 (clear vs. ambiguous) + ...
  indexDF <- data.frame(matrix(NA, nrow = length(tvec), ncol = (length(bloop)+3)))
  # bloop2 will become the column names of indexDF, so make it pretty:
  bloop2 <- bloop
  bloop2[str_which(bloop2, "; ")] <- substr(bloop2[str_which(bloop, "; ")], start = 1, str_length(bloop2)[str_which(bloop2, "; ")] - 2)
  bloop2 <- str_replace_all(bloop2, "; ", " or ")
  colnames(indexDF) <- c("Assigned", "Ambiguous", "Not.Assigned", bloop2)
  # popluate index df:
  indexDF$Assigned[1:length(i.unam)] <- i.unam
  indexDF$Ambiguous[1:length(i.am)] <- i.am
  indexDF$Not.Assigned[1:length(c(tvec.containa, which(tvec.na)))] <- sort(c(tvec.containa, which(tvec.na)))
  for (i in 1:length(bloop)) {
    boi <- which(tvec %in% bloop[i])
    indexDF[1:length(boi),bloop2[i]] <- boi
  }
  # create a plotting dataframe using length of non-NA columns in indexDF:
  plotDF.asv <- data.frame(trait.ass = colnames(indexDF),
                           nASV = apply(indexDF, MARGIN = 2, function(x) nrow(indexDF) - sum(is.na(x))), 
                           pASV = apply(indexDF, MARGIN = 2, function(x) (nrow(indexDF) - sum(is.na(x)))/nrow(indexDF)))
  # Plot mapping results in terms of Number (maybe change to proportion) of ASVs in each classification
  p1 <- ggplot(plotDF.asv[!plotDF.asv$trait.ass %in% c("Assigned", "Ambiguous"),], aes(x = trait.ass, y = pASV)) + 
    geom_bar(stat="identity", color = "black", position=position_dodge()) + 
    labs(x = "Trait Assignment", y = "Proportion of ASVs") + 
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
    ggtitle(paste0("Trait mapping results for ", trait.name))
  # if an otu-table was provided, compute the sample-wise proportion of seq reads mapped to each trait category
  if (!identical(otu.table, "none")) {
    if (any(rowSums(otu.table) > 1)){
      # compute relative read abundances if raw reads were supplied
      otu.table <- otu.table / rowSums(otu.table)
    }
    # return(list(otu.table, map.result, indexDF))
    raDF <- data.frame(matrix(NA, nrow = nrow(otu.table), ncol = ncol(indexDF))) # an nsample x ntrait.category matrix for storing relative abundances
    colnames(raDF) <- colnames(indexDF)
    # return(list(otu.table, indexDF, raDF))
    for (i in 1:ncol(indexDF)) {
      cidx <- indexDF[,i]
      cidx <- cidx[!is.na(cidx)]
      # cidx is the col indices of ASVs w/ the ith trait classification in indexDF
      if (length(cidx) > 1) {
        raDF[,i] <- rowSums(otu.table[,cidx])
      } else if (length(cidx) == 1) {
        raDF[,i] <- otu.table[,cidx]
      } 
    }
    # raDF now contains the summed relative read abundances of all ASVs within each trait classification (cols) for each sample (rows)
    
    # boxplot cumulative read abundances across all samples 
    plotDF2 <- melt(raDF)
    p2 <- ggplot(plotDF2[!plotDF2$variable %in% c("Assigned", "Ambiguous"),]) + 
      geom_boxplot(stat="boxplot", aes(x=variable, y=value), position=position_dodge()) + 
      labs(x = "Trait Assignment", y = "Cumulative Rel. Abundace Across Samples") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_text(size = 12, face="bold"),
            axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12, face="bold"),
            panel.background = element_rect(fill = "white",
                                            colour = "white",
                                            linetype = "solid"),
            panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                            colour = "white"),
            panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                            colour = "white"),
            axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
            legend.position = "none") + 
      stat_summary(fun = "mean", aes(x=variable, y=value), geom = "crossbar", colour = "red", width = 0.75)

    # histogram of rel. read abundances across samples that did not have definitive trait classifications
    p3 <- ggplot(raDF, aes(x = Ambiguous)) + 
      geom_histogram()  + 
      labs(x = "Rel. Read Abundance of ASVs w/ Unknown Trait Assignment", y = "N Samples") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_text(size = 12, face="bold"),
          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12, face="bold"),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "white"),
          axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"))
    
  } else {
    # all variables returned when otu.table is provided are returned as NULL
    raDF <- NULL
    plotDF2 <- NULL
    p2 <- NULL
    p3 <- NULL
  }
  return(list(p1, p2, p3, indexDF, plotDF.asv, raDF))
}


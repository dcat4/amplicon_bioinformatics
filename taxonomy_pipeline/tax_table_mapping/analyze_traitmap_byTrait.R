# UNDER Construction

# analyzes trait mapping outputs to determine the number of ASVs and distribution of read abundances/sample
# that can be mapped to each categorization within a user-specified trait. Optimized for categorical traits

# seems to work for number/proportion of ASVs (it runs for several traits I tested)
# need to check that ^that's accurate and add otu-table analysis

# map.result should be a df of output of traitmapper_Ramond
# trait.name is length 1 character vector with the name of trait (one colname of map.result) you want to analyze
# trait.options is a character vector with all possible assignments for the trait.name
# otu.table is a dataframe of sequence read counts for each ASV according to sample name 
# (^similar to phyloseq otu_table but should be a dataframe or matrix)

# outputs:
# 1. index array of ASVs with:
#   a. ambiguous trait assignments
#   b. unambiguous trait assignments
#   c. etc. each unique type of assignment for that trait (a breakdown of both am and unam assignments)
# 2. plots:
#   a. % of ASVs assigned to each trait, and ambiguous...
#   b. % of reads in each sample that is abiguous vs. unambiguous
#   c. % of reads in each sample assigned to each trait

analyze_traitmap_byTrait <- function(map.result, trait.name, trait.options, otu.table = "none", plotfilez = "none") {
  library("stringr")
  library("ggplot2")
  tvec <- map.result[,trait.name]
  tvec.na <- is.na(tvec) # where the name was mapped directly to NA, and nothing else
  tvec.containa <- str_which(tvec, "NA") # where the hit was NA in some instances but something else in others...
  nhitz <- str_count(tvec, ";")
  # nhitz will be 0 if there was 1 unique value for the trait
  # otherwise it will == the number of unique values of this trait
  
  # NOTE: need to add something that adds anything containing string "NA" in it's compiled traits to i.ambiguous (can't use which != 0 b/c of NA problems with which)
  
  # this works:
  i.am <- sort(c(which(tvec.na), which(nhitz > 0))) # ambiguous
  i.unam <- which(nhitz == 0) # unambiguous
  
  # return(list(tvec, tvec.na, nhitz, tvec.containa, i.unam, i.am))
  
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
  p1 <- ggplot(plotDF.asv[!plotDF.asv$trait.ass %in% c("Assigned", "Ambiguous"),], aes(x = trait.ass, y = nASV)) + 
    geom_bar(stat="identity", color = "black", position=position_dodge()) + 
    labs(x = "Trait Assignment", y = "Number of ASVs") + 
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
  
  
  
  return(list(p1,indexDF, plotDF.asv))
}


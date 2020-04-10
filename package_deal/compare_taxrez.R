# function originally written by Connie Liang
# modified and cleaned up by D. Catlett

# summarizes resolution of taxonomic assignments for any number of taxonomy tables of OTU/ASV and corresponding taxonomic assignments
# plots results as a stacked bar plot 

# function that returns a dataframe with proportion of ASVs of unassigned at each rank in each taxonomy table
# also makes and saves a grouped bar plot of proportion of ASVs unassigned at each rank in each taxonomy table
compare_taxrez <- function(..., pltfilez = "none", 
                           tablenames = c("bayes", "idtax"), 
                           ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species")){
  # package download/loading/updating:
  library("ggplot2") 
  library("reshape2")
  # a helper function for computing proportion unassigned
  nasum <- function(taxdf){
    notuz <- nrow(taxdf)
    x <- is.na(taxdf)
    ii <- colSums(x) / notuz
    return(ii)
  }
  
  x <- list(...)
  
  durp <- lapply(x, nasum)
  yaboi <- matrix(unlist(durp), nrow=length(tablenames), byrow=TRUE)
  colnames(yaboi) <- ranknamez
  rownames(yaboi) <- tablenames
  yaboi <- as.data.frame(t(yaboi))
  yaboi$rankz <- rownames(yaboi)
  
  yaboi <- melt(yaboi, id.vars = "rankz")
  p2 <- ggplot(yaboi, aes(fill = variable, x = rankz, y = value)) + 
    geom_bar(stat="identity", color = "black", position=position_dodge(width=0.8)) + 
    labs(x = "Taxonomic Rank", y = "Proportion Unassigned") + 
    scale_x_discrete(limits = ranknamez) + 
    scale_y_discrete(limits = seq(0.1,1.05,0.1), expand = c(0,0)) +
    coord_cartesian(ylim = c(0, 1)) + 
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
      scale_fill_discrete(name = "Taxonomy")
  if (pltfilez == "none"){
    # don't save anything
  } else {
    ggsave(filename = pltfile, plot = p2, device = "pdf")
  }
  return(list(yaboi, p2))
}

# under construction

# trying to do rank-wise comparisons for 3 tax tables

# almost there, just gotta add where all named and 2 agree + 1 don't

compare_byRank_3way <- function(table1, table2, table3,
                                pltfilez = c("prop_2wayplt.pdf", "abs_2wayplt.pdf"),
                                tablenames = c("bayes", "idtax"), 
                                ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species")) {
  notuz <- nrow(table1) # number of ASVs/OTUs/rows in each tax table
  # initialize a list for storing index vectors:
  allofit <- list()
  # you have to apply the below to each column, potentially in a loop
  for (i in 1:length(ranknamez)) {
    i.allna <- intersect(intersect(which(is.na(table1[,i])), which(is.na(table2[,i]))), which(is.na(table3[,i])))
    i.t1N.t23na <- intersect(intersect(which(!is.na(table1[,i])), which(is.na(table2[,i]))), which(is.na(table3[,i])))
    i.t2N.t13na <- intersect(intersect(which(is.na(table1[,i])), which(!is.na(table2[,i]))), which(is.na(table3[,i])))
    i.t3N.t12na <- intersect(intersect(which(is.na(table1[,i])), which(is.na(table2[,i]))), which(!is.na(table3[,i])))
    
    i.t12N.t3na <- intersect(intersect(which(!is.na(table1[,i])), which(!is.na(table2[,i]))), which(is.na(table3[,i])))
    i.t12N.t3na.a <- intersect(i.t12N.t3na, which(table1[,i] == table2[,i]))
    i.t12N.t3na.d <- intersect(i.t12N.t3na, which(table1[,i] != table2[,i]))
    
    i.t13N.t2na <- intersect(intersect(which(!is.na(table1[,i])), which(is.na(table2[,i]))), which(!is.na(table3[,i])))
    i.t13N.t2na.a <- intersect(i.t13N.t2na, which(table1[,i] == table3[,i]))
    i.t13N.t2na.d  <- intersect(i.t13N.t2na, which(table1[,i] != table3[,i]))
    
    i.t23N.t1na <- intersect(intersect(which(is.na(table1[,i])), which(!is.na(table2[,i]))), which(!is.na(table3[,i])))
    i.t23N.t1na.a <- intersect(i.t23N.t1na, which(table2[,i] == table3[,i]))
    i.t23N.t1na.d <- intersect(i.t23N.t1na, which(table2[,i] != table3[,i]))
    
    i.allN <- intersect(which(!is.na(table1[,i])), intersect(which(!is.na(table2[,i])), which(!is.na(table3[,i]))))
    i.allN.a <- intersect(i.allN, intersect(which(table1[,i] == table3[,i]), which(table2[,i] == table3[,i])))
    i.allN.d <- intersect(i.allN, intersect(which(table1[,i] != table3[,i]), intersect(which(table2[,i] != table3[,i]), which(table1[,i] != table2[,i]))))  # add another intersection where 1 != 2
    i.allN.12a <- intersect(i.allN, intersect(which(table1[,i] == table2[,i]), which(table2[,i] != table3[,i])))
    i.allN.13a <- intersect(i.allN, intersect(which(table1[,i] != table2[,i]), which(table1[,i] == table3[,i])))
    i.allN.23a <- intersect(i.allN, intersect(which(table1[,i] != table2[,i]), which(table2[,i] == table3[,i])))
    
    # prob merge these into a dataframe and make each entry in index list correspond to a rank...
    nr <- max(c(length(i.allna),length(i.t1N.t23na),length(i.t2N.t13na),length(i.t3N.t12na),
            length(i.t12N.t3na.a),length(i.t12N.t3na.d),length(i.t13N.t2na.a),length(i.t13N.t2na.d),length(i.t23N.t1na.a),length(i.t23N.t1na.d),
            length(i.allN.12a), length(i.allN.13a), length(i.allN.23a), length(i.allN.a), length(i.allN.d)))
    indexDF <-  data.frame(matrix(NA, nrow = nr, ncol = 15))
    colnames(indexDF) <- c("all.na", "t1named.t23NA", "t2named.t13NA", "t3named.t12NA", 
                           "t12named.same.t3NA", "t12named.diff.t3NA", "t13named.same.t2NA", "t13named.diff.t2NA", "t23named.same.t1NA", "t23named.diff.t1NA",
                           "allN.t12same.t3diff","allN.t13same.t2diff","allN.t23same.t1diff", "all.same.name", "all.diff.name")
    # all NA:
    if (length(i.allna) > 0){
      indexDF[1:length(i.allna), "all.na"] <- i.allna
    }
    # t1 named, 2/3 NA
    if (length(i.t1N.t23na) > 0){
      indexDF[1:length(i.t1N.t23na), "t1named.t23NA"] <- i.t1N.t23na
    }
    # t2 named, 1/3 NA
    if (length(i.t2N.t13na) > 0){
      indexDF[1:length(i.t2N.t13na), "t2named.t13NA"] <- i.t2N.t13na
    }
    # t3 named, 1/2 NA
    if (length(i.t3N.t12na) > 0){
      indexDF[1:length(i.t3N.t12na), "t3named.t12NA"] <- i.t3N.t12na
    }
    # t1/2 named agree t3 NA:
    if (length(i.t12N.t3na.a) > 0){
      indexDF[1:length(i.t12N.t3na.a), "t12named.same.t3NA"] <- i.t12N.t3na.a
    }
    if (length(i.t12N.t3na.d) > 0){
      indexDF[1:length(i.t12N.t3na.d), "t12named.diff.t3NA"] <- i.t12N.t3na.d
    }
    # t1/3 named agree t2 NA:
    if (length(i.t13N.t2na.a) > 0){
      indexDF[1:length(i.t13N.t2na.a), "t13named.same.t2NA"] <- i.t13N.t2na.a
    }
    if (length(i.t13N.t2na.d) > 0){
      indexDF[1:length(i.t13N.t2na.d), "t13named.diff.t2NA"] <- i.t13N.t2na.d
    }
    # t2/3 named agree t1 NA:
    if (length(i.t23N.t1na.a) > 0){
      indexDF[1:length(i.t23N.t1na.a), "t23named.same.t1NA"] <- i.t23N.t1na.a
    }
    if (length(i.t23N.t1na.d) > 0){
      indexDF[1:length(i.t23N.t1na.d), "t23named.diff.t1NA"] <- i.t23N.t1na.d
    }
    
    # all named in agreement:
    if (length(i.allN.a) > 0){
      indexDF[1:length(i.allN.a), "all.same.name"] <- i.allN.a
    }
    # all named in disagreement:
    if (length(i.allN.d) > 0){
      indexDF[1:length(i.allN.d), "all.diff.name"] <- i.allN.d
    }
    # t1/2 agree but t3 doesnt:
    if (length(i.allN.12a) > 0){
      indexDF[1:length(i.allN.12a), "allN.t12same.t3diff"] <- i.allN.12a
    }
    # t1/3 agree but t2 doesnt:
    if (length(i.allN.13a) > 0){
      indexDF[1:length(i.allN.13a), "allN.t13same.t2diff"] <- i.allN.13a
    }
    # t2/3 agree but t1 doesnt:
    if (length(i.allN.23a) > 0){
      indexDF[1:length(i.allN.23a), "allN.t23same.t1diff"] <- i.allN.23a
    }
    
    # below is copy/pasted from 2-way:
    pp <- apply(indexDF, MARGIN = 2, function(x) length(which(!is.na(x)))) # one rank's worth of numbers you can plot
    if (i == 1) {
      plotDF <- data.frame(comp = names(pp), count = pp, rank = rep(ranknamez[i], times = length(pp)))
    } else {
      plotDF <- rbind(plotDF, data.frame(comp = names(pp), count = pp, rank = rep(ranknamez[i], times = length(pp))))
    }
    # return(list(indexDF, pp, table1, table2))
    allofit[[i]] <- indexDF
  }
  # outside the loop now
  # move to plotting:
  plotDF[,"count"] <- plotDF[,"count"] / notuz
  
  # plotting with ggplot2:
  library("ggplot2")
  # proportional number of ASVs
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
    scale_fill_discrete(name = "Comparison", labels = c(all.na = "All 3 NA", 
                                                        t1named.t23NA = paste0(tablenames[1], " named,\n", tablenames[2], ", ", tablenames[3], " NA"), 
                                                        t2named.t13NA = paste0(tablenames[2], " named,\n", tablenames[1],", ", tablenames[3], " NA"), 
                                                        t3named.t12NA = paste0(tablenames[3], " named,\n", tablenames[1],", ", tablenames[2], " NA"), 
                                                        t12named.same.t3NA = paste0(tablenames[1],", ",tablenames[2], " same name,\n", tablenames[3], " NA"),
                                                        t13named.same.t2NA = paste0(tablenames[1],", ",tablenames[3], " same name,\n", tablenames[2], " NA"),
                                                        t23named.same.t1NA = paste0(tablenames[2],", ",tablenames[3], " same name,\n", tablenames[1], " NA"),
                                                        t12named.diff.t3NA = paste0(tablenames[1],", ",tablenames[2], " diff. name,\n", tablenames[3], " NA"),
                                                        t13named.diff.t2NA = paste0(tablenames[1],", ",tablenames[3], " diff. name,\n", tablenames[2], " NA"),
                                                        t23named.diff.t1NA = paste0(tablenames[2],", ",tablenames[3], " diff. name,\n", tablenames[1], " NA"),
                                                        allN.t12same.t3diff = paste0(tablenames[1],", ",tablenames[2], " same name,\n", tablenames[3], " diff. name"),
                                                        allN.t13same.t2diff = paste0(tablenames[1],", ",tablenames[3], " same name,\n", tablenames[2], " diff. name"),
                                                        allN.t23same.t1diff = paste0(tablenames[2],", ",tablenames[3], " same name,\n", tablenames[1], " diff. name"),
                                                        all.same.name = "All same name", all.diff.name = "All diff. name")) +
    ggtitle(paste0(tablenames[1], ", ",  tablenames[2], ", ",  tablenames[3]))
  
  # absolute number of ASVs:
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
    scale_x_discrete(labels = c(all.na = "All 3 NA", 
                                t1named.t23NA = paste0(tablenames[1], " named,\n", tablenames[2], ", ", tablenames[3], " NA"), 
                                t2named.t13NA = paste0(tablenames[2], " named,\n", tablenames[1],", ", tablenames[3], " NA"), 
                                t3named.t12NA = paste0(tablenames[3], " named,\n", tablenames[1],", ", tablenames[2], " NA"), 
                                t12named.same.t3NA = paste0(tablenames[1],", ",tablenames[2], " same name,\n", tablenames[3], " NA"),
                                t13named.same.t2NA = paste0(tablenames[1],", ",tablenames[3], " same name,\n", tablenames[2], " NA"),
                                t23named.same.t1NA = paste0(tablenames[2],", ",tablenames[3], " same name,\n", tablenames[1], " NA"),
                                t12named.diff.t3NA = paste0(tablenames[1],", ",tablenames[2], " diff. name,\n", tablenames[3], " NA"),
                                t13named.diff.t2NA = paste0(tablenames[1],", ",tablenames[3], " diff. name,\n", tablenames[2], " NA"),
                                t23named.diff.t1NA = paste0(tablenames[2],", ",tablenames[3], " diff. name,\n", tablenames[1], " NA"),
                                allN.t12same.t3diff = paste0(tablenames[1],", ",tablenames[2], " same name,\n", tablenames[3], " diff. name"),
                                allN.t13same.t2diff = paste0(tablenames[1],", ",tablenames[3], " same name,\n", tablenames[2], " diff. name"),
                                allN.t23same.t1diff = paste0(tablenames[2],", ",tablenames[3], " same name,\n", tablenames[1], " diff. name"),
                                all.same.name = "All same name", all.diff.name = "All diff. name"))
  ggtitle(paste0(tablenames[1], ", ",  tablenames[2], ", ",  tablenames[3]))
  
  if (length(grep(pltfilez, "none")) == 1){
    # don't save anything
  } else {
    ggsave(filename = pltfile, plot = p1, device = "pdf")
    ggsave(filename = pltfile, plot = p2, device = "pdf")
  }
  return(list(indexDF, plotDF, p1, p2))
}
  
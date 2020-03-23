# mines the pr2 database R package for all unique taxonomic paths. Saves output as a .csv file.

pr2_tax_mine <- function(out.file = "pr2_all_tax.csv") {
  library("pr2database")
  data("pr2")
  r <- c("kingdom", "supergroup", "division", "class", "order", "family", "genus", "species")
  tt <- pr2[, which(colnames(pr2) %in% r)]
  ty <- unique(tt, MARGIN = 1)
  if (length(grep(out.file, "none")) == 0) {
    write.csv(ty, out.file)
    return(ty)
  } else {
    return(ty)
  }
}


# this is a function that takes a consensus/ensemble taxonomy and the taxonomy tables it was built from and does a QC
# for frankenstein assignments. Hoping I can just write it in as a check at the end of all the consensus taxonomy algorithms

# I think it relies on dplyr's 
check4frankenstein <- function(tbl.list, ranknames = c("kingdom","supergroup","division","class","order","family","genus","species")) {

  library("dplyr")
  # pull out ensemble and the output of each individual algorithm:
  ee <- tbl.list[["ensemble"]]
  tt <- tbl.list[names(tbl.list) != "ensemble"]
  tt <- bind_rows(tt)
  
  # ee <- ee[, colnames(ee) %in% ranknames]
  # tt <- tt[, colnames(tt) %in% ranknames]
  # 
  # extract unique taxonomic paths from the ensemble and all individual tables:
  uee <- unique(ee, MARGIN = 1)
  tl <- unique(tt, MARGIN = 1)
  
  uee <- ee
  tl <- tt
  nchex <- nrow(uee)
  for (row in 1:nchex) {
    checker <- intersect(uee[row , ], tl)
    if (nrow(checker) == 0){
      ee <- uee[row , ]
      tmp <- tl
      tmp <- tmp[, 1:max(which(!is.na(ee)))]
      ee <- ee[, 1:max(which(!is.na(ee)))]
      check2 <- intersect(ee, tmp)
      if (nrow(check2) == 0){
        # this is a frankensteined assignment. throw an error and break
        stop("Frankenstein assignments detected. You need to modify your ensemble parameters and try again. \n Check Row ", toString(row))
      }
    }
  }
  
  # if you make it here, there are no frankenstein assignments. print that and you're done
  message("no Frankenstein assignments detected")
}



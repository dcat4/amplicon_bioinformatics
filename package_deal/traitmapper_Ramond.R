# function to map taxonomy table to the trait database published by Ramond et al 2019. 

# it becomes really slow and memory intensive if you compile all trait hitz when you're mapping
# broad names like "Eukaryota". I built in some work-arounds for this:
# 1. specify names you don't want to map w/ dont.map
# 2. it doesn't compile entries for taxonomic columns included in map2.taxnames, just uses NA

# taxin is the taxtable you want to map
# map2 is the trait database as a dataframe
# map2.taxnames is the column names of taxonomy entries in the trait database (the names you're mapping taxin names onto)
# NOTE: values in these columns are not compiled in the final mapped output array -- this greatly reduces memory usage and run time
# NOTE (cont): you can change this if it becomes necessary to look at all the mapped taxonomies...
# dont.map is a character vector of names to skip mapping -- use it for things like "Eukaryota" that will have 
# no meaningful trait annotations. Default is all unique Kingdom + Supergroup names from pr2 v.4.12.0
# filezout, if not "none", saves the mapped output and the names not mapped to RDS's or csv's

traitmapper_Ramond <- function(taxin, map2, 
                               map2.taxnames = c("Lineage1","Lineage2","Lineage3","Lineage4",
                                                 "Lineage5","Lineage6","Lineage7","Lineage8",
                                                 "Lineage9","Lineage10","Lineage11","Lineage12",
                                                 "Fam","Taxogroup","Taxo1","Last"),
                                dont.map = c("Eukaryota", "Archaea", "Bacteria", 
                                              "Alveolata","Opisthokonta","Archaeplastida","Excavata","Rhizaria","Stramenopiles",
                                              "Hacrobia","Amoebozoa","Apusozoa","Eukaryota_X","Protalveolata","Terrabacteria"),
                                filezout = "none") {
  bloop <- data.frame(matrix(NA, nrow = nrow(taxin), ncol = ncol(map2)+2))
  mapout <- cbind(taxin, bloop) # includes your taxonomy array appended to the mapping results (NA's right now)
  nomap <- c() # append names that don't hit the trait db here
  nranks <- ncol(taxin)
  counter <- nranks # for moving backwards down the tax table in mapping
  for (i in 1:nranks) {
    nn <- unique(taxin[,counter]) # unique names at this rank that you're going to map
    nn[nn %in% dont.map] <- NA # NA out names that you don't need to map 
    nn <- unique(nn) # unique names at this rank that you're going to map
    for (j in 1:length(nn)) {
      if (isFALSE(is.na(nn[j]))) {
        # you need to try to map this name...
        bloop <- apply(map2[,map2.taxnames], MARGIN = c(1,2), function(x) nn[j] == x)
        matchr <- which(rowSums(bloop, na.rm = TRUE) > 0) # rows that contain the name being mapped
        if (length(matchr) == 0) {
          # no hits - append name to nomap:
          nomap <- append(nomap, nn[j])
        } else if (length(matchr) == 1) {
          # there is one match. find rows of your mapout (combined taxonomy and trait db)
          # that have this name and pop the designated traits into the mapout array
          
          # find where mapout trait-mapped is na (hasn't been mapped yet) and where it == this name
          bloop2 <- intersect(which(is.na(mapout[,nranks+1])), which(mapout[,counter] == nn[j]))
          # pop it in to the mapped array:
          mapout[bloop2, (nranks+1):ncol(mapout)] <- as.data.frame(cbind(matrix(c(nn[j], length(matchr)), nrow=(length(bloop2)), ncol=2,byrow=TRUE),
                                                         matrix(map2[matchr,], nrow = (length(bloop2)), ncol = ncol(map2), byrow = TRUE)))
        } else if (length(matchr) > 1) {
          # there were multiple hits. find them & compile unique traits into one entry
          hitz <- map2[matchr,]
          mapme <- vector(mode = "character", length = ncol(hitz))
          for (k in 1:ncol(hitz)) {
            cn <- colnames(hitz)[k]
            uu <- unique(hitz[,k])
            if (length(uu) > 1 && !(cn %in% map2.taxnames)) {
              u2 <- paste0(uu[1],"; ")
              for (l in 2:length(uu)) {
                u2 <- paste0(u2, uu[l], "; ")
              }
              mapme[k] <- u2
            } else if (length(uu) > 1 && (cn %in% map2.taxnames)) {
              mapme[k] <- NA
            } else {
              mapme[k] <- uu
            }
          }
          # outside of the hitz/ k for-loop
          # within the matchr > 1 if statement
          # HERE you should find matches to nn[j] in mapout, and fill the trait mapping with mapme
          # find where mapout trait-mapped is na (hasn't been mapped yet) and where it == this name
          bloop2 <- intersect(which(is.na(mapout[,nranks+1])), which(mapout[,counter] == nn[j]))
          # pop it in to the mapped array:
          mapout[bloop2, (nranks+1):ncol(mapout)] <- as.data.frame(cbind(matrix(c(nn[j], length(matchr)), nrow=(length(bloop2)), ncol=2,byrow=TRUE),
                                                                         matrix(mapme, nrow = (length(bloop2)), ncol = ncol(map2), byrow = TRUE)),
                                                                   stringsAsFactors = FALSE)
        }
      }
    }
    # just outside the j for-loop
    counter <- counter - 1
  }
  # outside the i for-loop
  colnames(mapout) <- c(colnames(taxin),"name.mapped","n.hitz",colnames(map2))
  mapout <- as.data.frame(apply(mapout, MARGIN = 2, FUN = as.character), stringsAsFactors = FALSE)
  if (filezout == "none"){
    # don't save anything
    return(list(mapout, nomap))
  } else {
    write.csv(mapout, file = filezout[1])
    write.csv(nomap, file = filezout[2])
    return(list(mapout, nomap))
  }
}


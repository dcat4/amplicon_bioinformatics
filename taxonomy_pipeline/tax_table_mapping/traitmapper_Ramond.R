# UNDER CONSTRUCTION

# function to map taxonomy table to the trait database published by Ramond et al 2019. 

# should return a output file with the mapped results, and a second array of names that 
# were unable to be mapped

# got a decent start but going to need to see the datasets again to get it polished

# taxin is the taxtable you want to map
# map2 is the trait database as a dataframe
# map2.taxnames is the column names of taxonomy entries in the trait database (the names you're mapping taxin names onto)
# filezout, if not "none", saves the mapped output and the names not mapped to RDS's or csv's

traitmapper <- function(taxin, map2, 
                        map2.taxnames = c("Lineage1","Lineage2","Lineage3","Lineage4",
                                          "Lineage5","Lineage6","Lineage7","Lineage8",
                                          "Lineage9","Fam","Taxogroup","Taxo1","Last"),
                        filezout = "none") {
  bloop <- data.frame(matrix(NA, nrow = nrow(taxin), ncol = ncol(map2)+2))
  mapout <- cbind(taxin, bloop) # includes your taxonomy array appended to the mapping results (NA's right now)
  nomap <- c() # append names that don't hit the trait db here
  nranks <- ncol(taxin)
  
  counter <- nranks # for moving backwards down the tax table in mapping
  for (i in 1:nranks) {
    nn <- unique(taxin[,counter]) # unique names at this rank that you're going to map
    
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
          
          # DO THAT^^ HERE!!
        } else if (length(matchr) > 1) {
          # there were multiple hits. find them & compile unique traits into one entry
          hitz <- map2[matchr,]
          mapme <- vector(mode = "character", length = ncol(hitz))
          for (k in 1:ncol(hitz)) {
            uu <- unique(hitz[,k])
            if (length(uu) > 1) {
              u2 <- paste0(uu[1],"; ")
              for (l in 2:length(uu)) {
                u2 <- paste0(u2, uu[l], "; ")
              }
              mapme[k,] <- u2
            } else {
              mapme[k,] <- uu
            }
          }
          # outside of the hitz/ k for-loop
          # within the matchr > 1 if statement
          # HERE you should find matches to nn[j] in mapout, and fill the trait mapping with mapme
          # 
        }
        
      }
      # this is outside the if isFALSE(isna(nn)) bit, which means here nn[j] was NA
    }
    # just outside the j for-loop
    counter <- counter - 1
  }
  # outside the i for-loop
  
  return(list(mapout, nomap))
}
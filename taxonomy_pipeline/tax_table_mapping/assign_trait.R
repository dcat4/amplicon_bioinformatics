assign_trait <- function(taxin, trait.map, traits=c("SizeMin", "SizeMax", "Cover", "Shape", "Spicule", "Symmetry", 
                                                     "Polarity", "Colony", "Motility", "Chloroplast", "Plast_Origin", 
                                                     "Ingestion", "Behaviour", "Mutualistic_Host", "Symbiontic", 
                                                     "Symbiont_Location", "Host_Specialisation", "Symbiont_Specialisation", 
                                                     "Prey_Specialisation", "Mucilage", "Chemical_Signal", "Nutrient_Afinity", 
                                                     "Oxygen_Tolerance", "Salinity", "Temperature", "Depth", "Toxygenity", 
                                                     "Benthic_Phase", "Longevity", "Cyst_Spore", "Ploidy", "Genome_Size"), 
                         num.traits=1, separate=TRUE) {
  
  finializeTrait <- function(x, n) {
    # if there is a semicolon, resolve this
    if (grepl(";", x, fixed = TRUE)) {
      # check if NA exists in the string, if it does, then you're done
      if (str_detect(x, "NA")) {
        return(NA)
      } else {
        # no need to split, just find the nth occurence of the semicolon and substring before that
        idx <- gregexpr(";", x)[[1]]
        if (n > length(idx)) {
          return(x)
        } else {
          return(substr(x, 1, idx[n] - 1))
        }
      }
    } else { # if not, just return whatever since it's a single item
      if (!is.na(x) && x == "NA") {
        return (NA)
      } else {
        return(x)
      }
    }
  }
  
  trait.sub <- trait.map[, c("svN", "ASV", traits)]
  mapped <- merge(taxin[, c("svN", "ASV")], trait.sub, by=c("svN", "ASV"))
  
  for (row in 1:nrow(mapped)) {
    for (col in traits) {
      mapped[row, col] <- finializeTrait(mapped[row, col], num.traits)
    }
  }
  
  return(mapped)
}
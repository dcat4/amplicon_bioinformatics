assign_trait <- function(taxin, trait.map, traits=c("SizeMin", "SizeMax", "Cover", "Shape", "Spicule", "Symmetry", 
                                                     "Polarity", "Colony", "Motility", "Chloroplast", "Plast_Origin", 
                                                     "Ingestion", "Behaviour", "Mutualistic_Host", "Symbiontic", 
                                                     "Symbiont_Location", "Host_Specialisation", "Symbiont_Specialisation", 
                                                     "Prey_Specialisation", "Mucilage", "Chemical_Signal", "Nutrient_Afinity", 
                                                     "Oxygen_Tolerance", "Salinity", "Temperature", "Depth", "Toxygenity", 
                                                     "Benthic_Phase", "Longevity", "Cyst_Spore", "Ploidy", "Genome_Size"), 
                         binary=FALSE, want=rep(NA, length(traits)), separate=TRUE) {
  
  searchTraits <- function(x, wanted) {
    # if there is a semicolon, resolve this
    if (grepl(";", x, fixed = TRUE)) {
      return (NA)
    } else { # if not, just return whatever since it's a single item
      if (is.na(x)) {
        return (NA)
      }
      else if (x == "NA") {
        return (NA)
      }
      else if (x == wanted) {
        return (TRUE)
      } else {
        return(FALSE)
      }
    }
  }
  
  finializeTrait <- function(x) {
    # if there is a semicolon, resolve this
    if (grepl(";", x, fixed = TRUE)) {
      return (NA) # if there are multiple possibilities, then can't be resolved with NA
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
  
  num.traits <- length(traits)
  n.rows <- nrow(mapped)
  
  for (row in 1:n.rows) {
    for (col in 1:num.traits) {
      if (binary) {
        mapped[row, traits[col]] <- searchTrait(mapped[row, traits[col]], want[col])
      } else {
        mapped[row, traits[col]] <- finializeTrait(mapped[row, traits[col]])
      }
    }
  }
  
  return(mapped)
}
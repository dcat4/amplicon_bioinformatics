consensus_tax_bestRez <- function(..., tablenames = c("bayes", "idtax"), ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"),
                                  tiebreakz = "LCAlike") {
  x <- list(...)
  # majority rule system
  # use the assignment that is in agreement in the most taxonomy tables
  # row and rank wise basis majority rule system if the majority of the input taxonomy
  
  # try to avoid for loops (lapply?)
  
  # assume homoegenous ranking convention -> each input follow the same ranking
  
  # have an arbitrary number of input of taxonomy tables where they are all formatted the same way
  # row 1 corresponds to row 1 of another table 
  
  # assuming that you have three input tables...
  # for each table in table
  # row 1 would be ASV1 and would have a series of assignments 
  # for each of those name
  # if this name is shared by two or more of the tables, retain the name (majority)
  
  # for each row, start at the highest column first and iterate backwards to first column
  # check for equivalency -> if equal, pick one row and assign to the output
  # if two of them are equal -> assign either 
  # if > 50% of the matching  -> assign it to the upstream name 
  # < 50% -> use tie breaking scheme in consensus -> copy/paste 
  # store those row indicies and loop through a tie breaking loop 
  
  # consider A = tax, B = NA, C = NA 
  # have another input if NA is considered as a NA
  # -> if that input is true, value is NA
  # have some input where you have one option as NA as standard name -> assign upstream names
  
  # start treating NA's as names
  # do it by row
  # store it for the tie breakers if anything downstream doesn't match or use NA
}
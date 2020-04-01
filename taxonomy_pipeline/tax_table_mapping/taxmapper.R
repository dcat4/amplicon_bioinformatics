# UNDER CONSTRUCTION
# mapping algorithm for mapping a taxonomy table onto a taxonomy by exact name-matching, regardless of rank

# taxin is a dataframe containing ASV sequences and a set of taxonomic assignments
# tax2map2 can either be a filename or a dataframe of the taxonomic nomenclature you'd like to "translate" taxin into (use pr2_all_tax.csv)
# synonym.file is a filename or dataframe (use tax_synonyms_FINAL.csv in this folder to build) containing synonyms to check if exact name matching doesn't work
# outfilez is a character vector specifying names of output csv files. if none, do not write outputs to files

# function should always return a dataframe and a character vector, and potentially a second dataframe
# the first dataframe should include all unique taxonomic assignments (unique rows of taxin, excluding ASV seqs) and their corresponding mapped assignments
# the character vector should include all names within taxin that were unable to be mapped
# the optional second data frame should have the ASV seqs supplied in taxin with the mapped taxonomic assignments (structure should be identical to taxin except for taxonomic nomenclature)

taxmapper <- function(taxin, tax2map2, 
                      synonym.file = "tax_synonyms_FINAL.csv", 
                      outfilez = "none") {
  
}

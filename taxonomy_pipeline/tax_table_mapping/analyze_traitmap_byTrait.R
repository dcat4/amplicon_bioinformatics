# UNDER Construction

# trying to write something to analyze mapped outputs

# map.result should be a matrix or df of output of traitmapper_Ramond
# trait.name is length 1 character vector with the name of trait (one colname of map.result) you want to analyze
# otu.table is... a dataframe of sequence read counts for each ASV according to sample name 
# (^similar to phyloseq otu_table but should be a dataframe or matrix)

# outputs:
# 1. index array of ASVs with:
#   a. ambiguous trait assignments
#   b. unambiguous trait assignments
#   c. etc. each unique type of assignment for that trait (a breakdown of unambiguous assignments)
# 2. plots:
#   a. % of ASVs assigned to each trait, and ambiguous...
#   b. % of reads in each sample that is abiguous vs. unambiguous
#   c. % of reads in each sample assigned to each trait

analyze_traitmap_byTrait <- function(map.result, trait.name, otu.table = "none", plotfilez = "none") {
  
  tvec <- map.result$trait.name
  tvec.na <- is.na(tvec)
  nhitz <- str_count(tvec, ";")
  # nhitz will be 0 if there was 1 unique value for the trait
  # otherwise it will == the number of unique values of this trait
  
  return(list(tvec,tvec.na,nhitz))
  # i.ambiguous <- sort(c(which(tvec.na), which()))
  
  return(list(p1,p2,indexDF))
}
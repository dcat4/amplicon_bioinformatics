# UNDER Construction

# trying to write something to analyze mapped outputs

# map.result should be a df of output of traitmapper_Ramond
# trait.name is length 1 character vector with the name of trait (one colname of map.result) you want to analyze
# otu.table is a dataframe of sequence read counts for each ASV according to sample name 
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
  
  tvec <- map.result[,trait.name]
  tvec.na <- is.na(tvec) # where the name was mapped directly to NA, and nothing else
  tvec.containa <- str_which(tvec, "NA") # where the hit was NA in some instances but something else in others...
  nhitz <- str_count(tvec, ";")
  # nhitz will be 0 if there was 1 unique value for the trait
  # otherwise it will == the number of unique values of this trait
  
  # NOTE: need to add something that adds anything containing string "NA" in it's compiled traits to i.ambiguous (can't use which != 0 b/c of NA problems with which)
  
  i.ambiguous <- sort(c(which(tvec.na), which(nhitz > 0)))
  i.clear <- which(nhitz == 0)
  return(list(tvec,tvec.na,nhitz))

  return(list(p1,p2,indexDF))
}


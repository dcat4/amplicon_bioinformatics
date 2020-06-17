# directory where your "package_deal" folder is:
dd <- "~/Documents/R/amplicon_bioinformatics/package_deal/"

# source everything in there:

# trait and tax mapping:
source(paste0(dd,"traitmapper_Ramond.R"))
source(paste0(dd,"traitmapper_Ramond_fast.R"))
source(paste0(dd,"assign_trait.R"))
source(paste0(dd,"taxmapper.R"))
source(paste0(dd,"pr2_tax_miner.R"))
source(paste0(dd,"analyze_traitmap_byTrait.R"))

# helper fcns:
source(paste0(dd,"LCA2df.R"))
source(paste0(dd,"idtax2df_pr2.R"))
source(paste0(dd,"bayestax2df.R"))
source(paste0(dd,"idtax2df_silva.R"))
source(paste0(dd,"find_asvs_by_name.R"))

# ensemble tax algrithms:
source(paste0(dd,"consensus_tax_LCAlike.R"))
source(paste0(dd,"consensus_tax_bestRez.R"))
source(paste0(dd,"consensus_tax_mostCom.R"))

# tax table comparison functions
source(paste0(dd,"compare_byRank_2way.R"))
source(paste0(dd,"compare_byRank_3way.R"))
source(paste0(dd,"compare_assignments_2way.R"))
source(paste0(dd,"compare_assignments_3way.R"))
source(paste0(dd,"compare_by_tax_name.R"))
source(paste0(dd,"compare_taxrez.R"))





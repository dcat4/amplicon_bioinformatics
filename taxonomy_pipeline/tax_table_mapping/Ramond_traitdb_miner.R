# this script cleans up the Ramond trait dataset so it's a pretty dataframe that can be used for trait-mapping
# output saved as an RDS

# downloaded from https://www.seanoe.org/data/00405/51662/

# citations:
# 1. Ramond Pierre, Siano Raffaele, Sourisseau Marc (2018). Functional traits of marine protists. SEANOE. https://doi.org/10.17882/51662
# 2. Ramond Pierre, Sourisseau Marc, Simon Nathalie, Romac Sarah, Schmitt Sophie, Rigaut-Jalabert Fabienne, Henry Nicolas, De Vargas Colomban, Siano Raffaele (2019). Coupling between taxonomic and functional diversity in protistan coastal communities. Environmental Microbiology, 21(2), 730-749. Publisher's official version : https://doi.org/10.1111/1462-2920.14537, Open Access version : https://archimer.ifremer.fr/doc/00478/58964/

rm(list = ls())
setwd("~/Documents/R/amplicon_bioinformatics/taxonomy_pipeline/tax_table_mapping")
library("stringr")

rdb <- read.table("Ramond_protist_traits.csv", header = FALSE, sep = ";", stringsAsFactors = FALSE)

# unpack lineages:
ramond.vars <- rdb[1,] # variable names
rdb <- rdb[-1,] # remove header row
l <- rdb$V1
nranks <-max(str_count(l, "\\|") + 1)
derp <- vector(mode = "character", length(nranks))
for (i in 1:nranks) {
  derp[i] <- paste0("Lineage", toString(i))
}
ramond.vars <- ramond.vars[-1]
ramond.vars <- c(derp, ramond.vars)

l <- unlist(lapply(l, str_replace_all, "\\+", "_"))
newl <- data.frame(matrix(NA, nrow = length(l), ncol = nranks))
for (i in 1:length(l)) {
  x <- l[i] 
  x <- unlist(strsplit(x, "\\|"))
  newl[i,1:length(x)] <- x
}
rdb <- rdb[,-1] # remove lineage column

clean.rdb <- cbind(newl, rdb)
colnames(clean.rdb) <- ramond.vars

saveRDS(clean.rdb, file = "Ramond_traitdb_clean.RDS")

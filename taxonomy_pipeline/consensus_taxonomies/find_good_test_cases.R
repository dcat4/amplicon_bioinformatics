# this script compiles good test-cases for evaluating consensus_most_Com

rm(list=ls())
library(dplyr)

# setwd and read in your datasets:
setwd("~/Documents/R/amplicon_bioinformatics/tax_pipe_Mar20/")
dd <- "~/Documents/R/amplicon_bioinformatics/package_deal/"

# read in all our functions:
source(paste0(dd,"all_of_it.R"))

# read in  the data
bayes.pr2 <- readRDS("initial_tax_tabs/bayes_pr2_0boot_Mar20.rds")
idtax.pr2 <- readRDS("initial_tax_tabs/idtax_pr2_0boot_Mar20.rds")
lca.pr2 <- readRDS("initial_tax_tabs/LCA_pr2_rawdf_Mar20.rds")

library("DECIPHER")

rubber <- readDNAStringSet("blaster/allASVs4blast.fasta")
# convert idtax and bayes to dataframes - keep boot = 0  and return conf for now...
# note that for idtaxa, order of ASVs must be retrieved from dada2 seqtab - rubric represents that order
# and adds ASV names and sequences to output dataframe:
xx <- bayestax2df(bayes.pr2, boot = 60, rubric = rubber, return.conf = TRUE)
bayes.pr2 <- xx[[1]]
bayes.pr2.conf <- xx[[2]]

xx <- idtax2df_pr2(idtax.pr2, boot = 50, rubric = rubber, return.conf = TRUE)
idtax.pr2 <- xx[[1]]
idtax.pr2.conf <- xx[[2]]

# re-order your 6 tables by sorting according to ASV:
ii <- base::sort(bayes.pr2$ASV, index.return = TRUE)
bayes.pr2 <- bayes.pr2[ii$ix,]
bayes.pr2.conf <- bayes.pr2.conf[ii$ix,]
ii <- base::sort(idtax.pr2$ASV, index.return = TRUE)
idtax.pr2 <- idtax.pr2[ii$ix,]
idtax.pr2.conf <- idtax.pr2.conf[ii$ix,]

tt <- read.csv(file = "~/Documents/R/amplicon_bioinformatics/taxonomy_pipeline/tax_table_mapping/pr2_all_tax.csv", stringsAsFactors = FALSE)
tt <- tt[, -c(1)]
# map lca.pr2 to pr2 taxonomy:
lca.pr2 <- taxmapper(lca.pr2, tax2map2 = tt, 
                     ignore.format = TRUE, exceptions = c("Bacteria", "Archaea"),
                     synonym.file = "~/Documents/R/amplicon_bioinformatics/taxonomy_pipeline/tax_table_mapping/tax_synonyms_FINAL.csv", outfilez = "none")
lca.pr2 <- lca.pr2[[3]]
ii <- base::sort(lca.pr2$ASV, index.return = TRUE)
lca.pr2 <- lca.pr2[ii$ix,]

# homogenize colnames:
colnames(idtax.pr2) <- colnames(bayes.pr2)
colnames(lca.pr2) <- colnames(bayes.pr2)

## just stepping thru different scenarios and manually curating a test set 
# find equivalent rows across all 3:
x <- dplyr::intersect(dplyr::intersect(bayes.pr2, idtax.pr2), lca.pr2)
x <- x[!is.na(x$Division) , ]
all.agree <- c("sv20747", "sv20858")

x <- dplyr::setdiff(dplyr::intersect(bayes.pr2, idtax.pr2), lca.pr2)
x <- x[!is.na(x$Division) , ]
lca.dis <- c("sv14136", "sv17278", "sv20738")

lca.pr2[lca.pr2$svN %in% all.agree , ]
idtax.pr2[idtax.pr2$svN %in% all.agree , ]
bayes.pr2[bayes.pr2$svN %in% all.agree , ]

lca.pr2[lca.pr2$svN %in% lca.dis , ]
idtax.pr2[idtax.pr2$svN %in% lca.dis , ]
bayes.pr2[bayes.pr2$svN %in% lca.dis , ]

x <- dplyr::setdiff(dplyr::setdiff(bayes.pr2, idtax.pr2), lca.pr2)
x <- x[!is.na(x$Division) , ]
all.dis <- c("sv3579", "sv20673", "sv8382")

lca.pr2[lca.pr2$svN %in% all.dis , ]
idtax.pr2[idtax.pr2$svN %in% all.dis , ]
bayes.pr2[bayes.pr2$svN %in% all.dis , ]

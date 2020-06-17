# this script analyzes the composition of the mock v9 data set I put together, look at the number of ASVs from different taxonomic groups
# and how variable the resolution of the expected assignments are.

# got some pie charts but there's a whole bunch of taxonomic groups (which is a good thing) that makes them illegible
# maybe just compare to the total number of lineages included in pr2

rm(list=ls())
setwd("~/Documents/R/amplicon_bioinformatics/mock_analysis/")

library("stringr")
library("ggplot2")

et <- readRDS("mock_data/v9_mock_exptax_bothPrimers.rds")
et[et == "NaN"] <- NA # because this didn't happen automatically...

# count NA's:
nact <- apply(et, 1, FUN = function(x) sum(is.na(x)))

hist(nact)

ii <- 4
dd <- as.data.frame(table(et[,ii]))
ggplot(dd, aes(x="", y=Freq, fill=Var1)) + 
  geom_bar(width=1, stat="identity") +
  coord_polar("y") +
  ggtitle(paste0("Unique ",colnames(et)[ii])) +
  theme_void() +
  theme(plot.title = element_text(hjust=0.5), legend.title=element_blank())

ii <- 5
dd <- as.data.frame(table(et[,ii]))
ggplot(dd, aes(x="", y=Freq, fill=Var1)) + 
  geom_bar(width=1, stat="identity") +
  coord_polar("y") +
  ggtitle(paste0("Unique ",colnames(et)[ii])) +
  theme_void() +
  theme(plot.title = element_text(hjust=0.5), legend.title=element_blank())

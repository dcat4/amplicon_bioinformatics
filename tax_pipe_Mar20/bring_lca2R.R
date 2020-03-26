# uses my helper function to bring LCA csv outputs into R, format them for use with the rest of the pipeline, and save as RDS's. 

rm(list=ls())

library("DECIPHER")

setwd("~/Documents/R/amplicon_bioinformatics/tax_pipe_Mar20")

silva.csv <- read.csv("initial_tax_tabs/MEGAN_LCA_Silva_Mar20.csv", header = FALSE, stringsAsFactors = FALSE)

p1.csv <- read.csv("initial_tax_tabs/MEGAN_LCA_pr2_half1.csv", header = FALSE, stringsAsFactors = FALSE)
p2.csv <- read.csv("initial_tax_tabs/MEGAN_LCA_pr2_half2.csv", header = FALSE, stringsAsFactors = FALSE)

pr2.csv <- rbind(p1.csv, p2.csv)

rubber <- readDNAStringSet("blaster/allASVs4blast.fasta")

source("~/Documents/R/amplicon_bioinformatics/taxonomy_pipeline/helper_fcns/LCA2df.R")
pr2.df <- LCA2df(pr2.csv, rubber)
silva.df <- LCA2df(silva.csv, rubber)

# save outputs:
saveRDS(pr2.df, file = "initial_tax_tabs/LCA_pr2_rawdf_Mar20.rds")
saveRDS(silva.df, file = "initial_tax_tabs/LCA_silva_rawdf_Mar20.rds")


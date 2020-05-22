# this script reformats the expected taxonomy data so that you have a mock-sequence-table (just the ASVs with all 1 counts that you can input to the 
# taxonomic assignment algorithms), saves that as an RDS.

# it also creates a 'rubric' (ASV sequences plus svN's), saves that as a fasta file for Blast/LCA analysis

# last it saves the expected taxonomy table to an rds for convenience

rm(list=ls())
setwd("~/Documents/R/amplicon_bioinformatics/mock_analysis/mock_data")

library("DECIPHER")
library("Biostrings")
library("stringr")

xx <- read.csv(file = "pr2_v9amp_exptax_bothPrimers.csv", stringsAsFactors = FALSE)

# clean up and save an expected taxonomy table:
y <- seq(from = 1, to = nrow(xx), by = 1)
y <- sapply(y, function(x) paste0("sv",toString(x)))
xx <- cbind(y,xx, stringsAsFactors = FALSE)
colnames(xx) <- c("svN","asv","exp.kingdom", "exp.supergroup", "exp.division", "exp.class", "exp.order", "exp.family", "exp.genus", "exp.species")

saveRDS(xx, file = "v9_mock_exptax_bothPrimers.rds")

# create and save a fake seqtab (same format as what dada2 outputs)
mat1 <- matrix(data=NA, nrow = 3, ncol = nrow(xx))
mat1[is.na(mat1)] <- 10
rownames(mat1) <- c("s1", "s2", "s3")
colnames(mat1) <- xx$asv

saveRDS(mat1, file = "v9_mock_seqtab_bothPrimers.rds")

# create and save a fastafile of mock asvs to do blast/LCA
s <- xx$asv
names(s) <- xx$svN
ff <- DNAStringSet(s, use.names = TRUE)
writeXStringSet(ff, filepath = "v9_mock_asvs_bothPrimers.fasta", format = "fasta")
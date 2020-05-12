# this script mines pr2 and silva databases to assemble a large test data set of V9 and V4 sequences...

# just copy/pasted, all of below, need to tweak

rm(list=ls())
setwd("~/Documents/R/amplicon_bioinformatics/mock_analysis/")
source("~/Documents/R/amplicon_bioinformatics/package_deal/all_of_it.R")

library("pr2database")
library("DECIPHER")
library("stringr")

data("pr2")

# the fields are described here:
# https://pr2-database.org/documentation/pr2-fields/

# extract reference sequences (according to pr2 creator)
refz <- pr2[!is.na(pr2$reference_sequence),]
# this chex out:
# > unique(refz$reference_sequence)
# [1] 1
fastaFile <- DNAStringSet(readRNAStringSet("~/Documents/R/silva_nr_v138_train_set.fa"))
# ^this throws a warning but i read it into matlab w/ no issues and same results
seq_name = names(fastaFile)
sequence = paste(fastaFile)
silva <- data.frame(seq_name, sequence, stringsAsFactors = FALSE)

# clear out unnecessary stuff:
rm("pr2", "fastaFile", "seq_name", "sequence")

# pull sequences from refz, fix up names to match silva, and combine silva/pr2
# pr2 names should go 
refz <- subset(refz, select = c("pr2_accession","species","kingdom","supergroup","division","class","order","family","genus","sequence"))

df <- as.character(NA, ncol = 1, nrow = nrow(refz))
for (i in 1:nrow(refz)) {
  # xx <- paste0(refz$pr2_accession[i], " ", refz$kingdom[i], "; ", refz$supergroup[i], "; ", refz$division[i], "; ", refz$class[i], "; "
  #              , refz$order[i], "; ", refz$family[i], "; ", refz$genus[i], "; ", refz$species[i])
  xx <- paste0(refz$pr2_accession[i],"|", refz$kingdom[i],"|", refz$supergroup[i],"|", refz$division[i],"|", refz$class[i],"|",
               refz$order[i],"|", refz$family[i],"|", refz$genus[i], "|",refz$species[i])
  df[i] <- xx
}
# worked, hooray! now df is a character vector w/ pr2 accession and taxonomy. 

pr2 <- data.frame(df, stringsAsFactors = FALSE)
pr2$sequence <- refz$sequence
colnames(pr2) <- colnames(silva)

# remove sequences in silva that are found in pr2 to make sure your ref data is all unique (this favors pr2)
ii <- which(silva$sequence %in% pr2$sequence)
i2 <- which(pr2$sequence %in% silva$sequence)

# shockingly, there are no overlapping sequences between silva + pr2. there are a couple hundred duplicate sequences in each though.
# remove duplicate sequences from each:

# merge pr2 and silva
both <- base::rbind(pr2, silva, stringsAsFactors = FALSE)

# create DNA string set object from merged:

# use one of matchProbePair, matchPattern, or pairwiseAlignment to extract amplicons


feck
# combine the two sets:
# there are 52767 sequences removed when you do unique(both$sequence).
# according to J. Bowman, that's bad and you need to remove them...
# he also says u should be using 1 gene at a time...
# see here: https://www.polarmicrobes.org/phylogenetic-placement-re-re-visited/
# I guess we can leave them in...

# subset for developing code:
sub <- both[1:500,]

# save a Fasta with your reference data set:
ref.fasta = dataframe2fas(sub, file="TEST_pr2_silva_merged_refdata.fasta")
nn <- sub$seq_name

# msa:
# https://bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf
# clustal omega alignment:
msa.co <- msa("TEST_pr2_silva_merged_refdata.fasta", method = "ClustalOmega", type = "dna", order = "input")
# muscle alignment, but u can't preserve order and muscle apparently is better for AA's anyhow
# msa.mu <- msa("TEST_pr2_silva_merged_refdata.fasta", method = "Muscle", type = "dna", order = "input")

msa.co <- unmasked(msa.co)
names(msa.co) <- nn # order was preserved with with 'order' argument in msa call so u can plug names in

# save msa as fasta:
writeXStringSet(msa.co, file = "TEST2_refdata_msa.fasta")

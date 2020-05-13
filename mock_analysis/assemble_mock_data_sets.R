# this script mines pr2 and silva databases to assemble a large test data set of V9 and V4 sequences...

# working on isolating amplicons
# running into issues w/ our primer set and silva database - almost never hits any of the ref seqs

rm(list=ls())
setwd("~/Documents/R/amplicon_bioinformatics/mock_analysis/")
source("~/Documents/R/amplicon_bioinformatics/package_deal/all_of_it.R")

library("pr2database")
library("DECIPHER")
library("Biostrings")
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

silvaeuk <- silva[which(str_detect(silva$seq_name,"Eukaryota")),]
hist(str_length(silvaeuk$sequence))

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
pr2 <- pr2[!(duplicated(pr2$sequence) | duplicated(pr2$sequence, fromLast = TRUE)), ]
silva <- silva[!(duplicated(silva$sequence) | duplicated(silva$sequence, fromLast = TRUE)), ]


### I'm here in writing.... ------------------

# isolate v9, v4, and v4-5 regions
# adapted from here: https://github.com/pr2database/pr2-primers/blob/master/PR2%20Primers%20pr2_match.R
pp <- read.csv("primer_sets.csv", stringsAsFactors = FALSE)
ur <- unique(pp$X18s_region)

pr2.2 <- pr2$sequence
names(pr2.2) <- pr2$seq_name
pr2 <- DNAStringSet(pr2.2)

silva.2 <- silva$sequence
names(silva.2) <- silva$seq_name
silva <- DNAStringSet(silva.2)

mm <- 2 # max mismatches for finding primer hits

pr2.prime <- vector(mode = "list", length = length(ur)*2)
silva.prime <- vector(mode = "list", length = length(ur)*2)
counter <- 1
for (i in 1:length(ur)) {
  fwd <-  DNAString(pp$sequence[min(which(pp$X18s_region == ur[i]))])
  rev <-  DNAString(pp$sequence[max(which(pp$X18s_region == ur[i]))])
  rev <-  reverseComplement(rev)

  fwd.pos.pr2 <- vmatchPattern(fwd, pr2, max.mismatch=mm, min.mismatch=0, with.indels=FALSE, fixed=FALSE, algorithm="auto")
  fwd.pos.silva <- vmatchPattern(fwd, silva, max.mismatch=mm, min.mismatch=0, with.indels=FALSE, fixed=FALSE, algorithm="auto")
  rev.pos.pr2 <- vmatchPattern(rev, pr2, max.mismatch=mm, min.mismatch=0, with.indels=FALSE, fixed=FALSE, algorithm="auto")
  rev.pos.silva <- vmatchPattern(rev, silva, max.mismatch=mm, min.mismatch=0, with.indels=FALSE, fixed=FALSE, algorithm="auto")
  
  # unlist the MIndex position
  fwd.pos.pr2 <- unlist(fwd.pos.pr2)
  fwd.pos.silva <- unlist(fwd.pos.silva)
  rev.pos.pr2 <- unlist(rev.pos.pr2)
  rev.pos.silva <- unlist(rev.pos.silva)

  # Create data frames to store primer match data
  # pr2:
  fhits.pr2 <- data.frame(namer = fwd.pos.pr2@NAMES, 
                         fpos.strt = fwd.pos.pr2@start, 
                         fpos.end = fwd.pos.pr2@start + (length(fwd)-1),
                         stringsAsFactors = FALSE)
  rhits.pr2 <- data.frame(namer = rev.pos.pr2@NAMES, 
                         rpos.strt = rev.pos.pr2@start,
                         rpos.end = rev.pos.pr2@start + (length(rev)-1),
                         stringsAsFactors = FALSE)
  # silva:
  fhits.silva <- data.frame(namer = fwd.pos.silva@NAMES, 
                           fpos.strt = fwd.pos.silva@start, 
                           fpos.end = fwd.pos.silva@start + (length(fwd)-1),
                           stringsAsFactors = FALSE)
  rhits.silva <- data.frame(namer = rev.pos.silva@NAMES, 
                           rpos.strt = rev.pos.silva@start,
                           rpos.end = rev.pos.silva@start + (length(rev)-1),
                           stringsAsFactors = FALSE)
  # store primer hit data in a list:
  pr2.prime[[counter]] <- fhits.pr2
  names(pr2.prime)[counter] <- paste0("fwd.",ur[i])
  silva.prime[[counter]] <- fhits.silva
  names(silva.prime)[counter] <- paste0("fwd.",ur[i])
  counter <- counter + 1
  pr2.prime[[counter]] <- rhits.pr2
  silva.prime[[counter]] <- rhits.silva
  counter <- counter + 1
}
# isolate 5000 sequences


# create DNA string set object from merged:


# use one of matchProbePair, matchPattern, or pairwiseAlignment to extract amplicons

feck

# merge pr2 and silva
both <- base::rbind(pr2, silva, stringsAsFactors = FALSE)

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

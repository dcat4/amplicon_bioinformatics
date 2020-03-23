# this script creates a trainingSet compatible with idtaxa algorithm for the pr2 database (output is saved)
# and then runs the idtaxa algorithm against the pr2 training set and saves the output

# bootstrapping threshold is set to 0 -- bootstrapping confidence values are returned in the saved output
# so that you can retroactively determine thresholds.

# creating the trainingSet for pr2:
dna <- readDNAStringSet("databases/pr2_version_4.12.0_18S_dada2.fasta") #change the path to your local server
tt <- names(dna)
taxonomy <- paste("Root", tt, sep=";")

trainingSetpr2 <- LearnTaxa(dna, taxonomy)
saveRDS(trainingSetpr2,"databases/trainSetpr2_Mar20.rds")

#idtaxa-ing
#assuming you have bioconductors and decipher packages

dna <- DNAStringSet(getSequences(seqtab))  #extracting the dna strings

ids <- IdTaxa(dna, trainingSetpr2, strand="both", threshold = 0) #the trainingSet here is what I have been working on, this trainingSet was specific to Silva, and we need the one for pr2
saveRDS(ids,"initial_tax_tabs/idtax_pr2_0boot_Mar20.rds")

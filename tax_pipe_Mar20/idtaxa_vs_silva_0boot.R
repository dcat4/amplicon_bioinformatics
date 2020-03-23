# idtaxa for silva
# bootstrapping threshold set to 0

dna <- DNAStringSet(getSequences(seqtab))

load("databases/SILVA_SSU_r138_2019.RData")

ids <- IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=FALSE, threshold = 0)

saveRDS(ids,"initial_tax_tabs/idtax_silva_0boot_Mar20.rds")


### this script does taxonomic assignments with the Bayesian classifier with no boot-strapping threshold

## Note: Assigning taxonomy takes a reeeeeeeeally long time.
tax.silva <- assignTaxonomy(seqtab, "databases/silva_nr_v138_train_set.fa", multithread=FALSE, tryRC=T, outputBootstraps=TRUE, minBoot=0)

## Save taxonomic assignments:
saveRDS(tax.silva, "initial_tax_tabs/bayes_silva_0boot_Mar20.rds") # CHANGE ME ...

tax.pr2 <- assignTaxonomy(seqtab, "databases/pr2_version_4.12.0_18S_dada2.fasta", multithread=TRUE, tryRC=T, outputBootstraps=TRUE, minBoot=0, taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"))

## Save taxonomic assignments:
saveRDS(tax.pr2, "initial_tax_tabs/bayes_pr2_0boot_Mar20.rds")






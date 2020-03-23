### chunks out the fasta file into halves b/c pr2 blast output is >10GB

rm(list=ls())

# set working directory
.cran_packages <- c("gridExtra", "knitr", "data.table", "ggplot2")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
    install.packages(.cran_packages[!.inst], repos='http://cran.us.r-project.org')
}

sapply(c(.cran_packages), require, character.only = TRUE)
library("DECIPHER")

setwd("~/Documents/R/desktop_ampData_processing/connie_taxonomy_stuff_Mar2020/18sV9_amplicon_sequencing/tax_pipe_Mar20/blaster")

seqer <- readDNAStringSet("allASVs4blast.fasta")

half1 <- seqer[1:12000,]
half2 <- seqer[12001:length(seqer),]

writeXStringSet(half1, filepath = "half1ASVs4blast.fasta", format = "fasta")
writeXStringSet(half2, filepath = "half2ASVs4blast.fasta", format = "fasta")
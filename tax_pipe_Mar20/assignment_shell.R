### this script does taxonomic assignments with the Bayesian classifier with no boot-strapping threshold

rm(list=ls())

# set working directory
.cran_packages <- c("gridExtra", "knitr", "data.table", "ggplot2")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
    install.packages(.cran_packages[!.inst], repos='http://cran.us.r-project.org')
}

sapply(c(.cran_packages), require, character.only = TRUE)
library("dada2")
library("DECIPHER")

setwd("/home/dcatlett/amp_data_processing/tax_pipe_Mar20")

seqtab <- readRDS("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/all18s_Mar20/chimera_removal/seqtab.rds")

# assignment scripts
source("idtaxa_vs_silva_0boot.R")
source("idtaxa_vs_pr2_0boot.R")
source("assign_bayesTax_0boot.R")
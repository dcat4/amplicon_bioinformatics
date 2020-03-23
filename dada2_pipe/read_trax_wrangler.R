# this script grabs the read count data from each step of dada 2 and merges it into a single .csv file
# this allows for easy chex of the number of reads passing each step of the pipeline.

rm(list=ls())

# set working directory
setwd("/home/dcatlett/amp_data_processing/dada2_ML_Mar20")
.cran_packages <- c("gridExtra", "knitr", "data.table", "ggplot2")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
    install.packages(.cran_packages[!.inst], repos='http://cran.us.r-project.org')
}

sapply(c(.cran_packages), require, character.only = TRUE)
library("dada2")
##### RE-establish your paths for chimera removal and onwards:

# write the name of your library here:
Lib1 <- "MockComm"
Lib2 <- "PBEuk01" 
Lib3 <- "PBEuk02" 
Lib4 <- "PBEuk03" 
Lib5 <- "PBEuk04" 
Lib6 <- "PBEuk05" 
Lib7 <- "TM01" 
Lib8 <- "TM02"
Lib9 <- "PBEuk06"
Lib10 <- "ACEuk01"

libnames <- c(Lib1, Lib2, Lib3, Lib4, Lib5, Lib6, Lib7, Lib8, Lib9, Lib10)
for (i in 1:length(libnames)) {
    filt <- read.csv(paste0("~/amp_data_processing/dada2_ML_Mar20/",libnames[i],"/processing_results/filtering/readsin_readsout.csv"), stringsAsFactors = FALSE)
    nami <- sapply(filt[,"X"], function(x) min(gregexpr("_",x)[[1]]) - 1)
    snam <- sapply(filt[,"X"], function(x) substr(x, 1, min(gregexpr("_",x)[[1]]) - 1))
    rownames(filt) <- snam

    samin <- read.csv(paste0("~/amp_data_processing/dada2_ML_Mar20/",libnames[i],"/sample_inference/readsMergedAndTabled.csv"), stringsAsFactors = FALSE)
    rownames(samin) <- samin$X
    this1 <- merge(filt, samin, by = "row.names")
    if (i == 1) {
        alld <- this1
    } else {
        alld <- rbind(alld, this1)
    }
    
}

rownames(alld) <- alld$Row.names
nochime <- read.csv("~/amp_data_processing/dada2_ML_Mar20/all18s_Mar20/chimera_removal/readsChimeraFiltered.csv", stringsAsFactors = FALSE)
rownames(nochime) <- nochime$X
alld <- merge(alld, nochime, by = "row.names")

rownames(alld) <- alld[,"Row.names"]
alld <- alld[,!(colnames(alld) %in% c("Row.names", "X.y", "X"))]
write.csv(alld, file = "readtrax_filt2chimremoval_Mar20.csv")

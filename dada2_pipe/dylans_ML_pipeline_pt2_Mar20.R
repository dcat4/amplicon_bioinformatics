### This version is for use with the ERI servers.

# This is my pipeline, adapted for use with multiple libraries. Basically just copy/pasting everything and numbering it...
# this is part 2 and picks up after inferring ASV's.
# this is how you clear your workspace:
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
Lib1<- "MockComm"
Lib2 <- "PBEuk01" 
Lib3 <- "PBEuk02" 
Lib4 <- "PBEuk03" 
Lib5 <- "PBEuk04" 
Lib6 <- "PBEuk05" 
Lib7 <- "TM01" 
Lib8 <- "TM02"
Lib9 <- "PBEuk06"
Lib10 <- "ACEuk01"

filtMethod <- ""
ISVmethod <- ""

pathout_ISV1 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib1,"/sample_inference/",filtMethod,ISVmethod)
pathout_ISV2 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib2,"/sample_inference/",filtMethod,ISVmethod)
pathout_ISV3 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib3,"/sample_inference/",filtMethod,ISVmethod)
pathout_ISV4 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib4,"/sample_inference/",filtMethod,ISVmethod)
pathout_ISV5 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib5,"/sample_inference/",filtMethod,ISVmethod)
pathout_ISV6 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib6,"/sample_inference/",filtMethod,ISVmethod)
pathout_ISV7 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib7,"/sample_inference/",filtMethod,ISVmethod)
pathout_ISV8 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib8,"/sample_inference/",filtMethod,ISVmethod)
pathout_ISV9 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib9,"/sample_inference/",filtMethod,ISVmethod)
pathout_ISV10 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib10,"/sample_inference/",filtMethod,ISVmethod)

###------------------------------------
## Step 1
## Combine libraries and remove chimeras with dada2_ML_remove_chimera_Nov2018_dc.R

projectName <- "all18s_Mar20"

# Where you want the output sequence table saved:
RCmethod <- "" # same as above, there's only 2 options here though really.
pathout_RC <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",projectName,"/chimera_removal/",filtMethod,ISVmethod)

source("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/dada2_ML_remove_chimera_Mar20.R")


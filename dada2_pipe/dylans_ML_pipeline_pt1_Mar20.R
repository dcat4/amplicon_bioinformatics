### This version is for use with the ERI servers.

# This is my pipeline, adapted for use with multiple libraries. Basically just copy/pasting everything and numbering it...

#### Outline of the dada2 amplicon sequence processing steps
### This version was put together by D. Catlett starting Mar2020

#### Note: Change directories below accordingly for each processing step and new library/parameter modification

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

##### Set up new directories to raw data and for storing results at each step of data processing:

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

# establish paths to raw miseq files
path1 <- paste0("/home/dcatlett/OCEANCOLOR/users/dcatlett/ampData_18sV9/Raw_data/",Lib1)
path2 <- paste0("/home/dcatlett/OCEANCOLOR/users/dcatlett/ampData_18sV9/Raw_data/",Lib2)
path3 <- paste0("/home/dcatlett/OCEANCOLOR/users/dcatlett/ampData_18sV9/Raw_data/",Lib3)
path4 <- paste0("/home/dcatlett/OCEANCOLOR/users/dcatlett/ampData_18sV9/Raw_data/",Lib4)
path5 <- paste0("/home/dcatlett/OCEANCOLOR/users/dcatlett/ampData_18sV9/Raw_data/",Lib5)
path6 <- paste0("/home/dcatlett/OCEANCOLOR/users/dcatlett/ampData_18sV9/Raw_data/",Lib6)
path7 <- paste0("/home/dcatlett/OCEANCOLOR/users/dcatlett/ampData_18sV9/Raw_data/",Lib7)
path8 <- paste0("/home/dcatlett/OCEANCOLOR/users/dcatlett/ampData_18sV9/Raw_data/",Lib8)
path9 <- paste0("/home/dcatlett/OCEANCOLOR/users/dcatlett/ampData_18sV9/Raw_data/",Lib9)
path10 <- paste0("/home/dcatlett/OCEANCOLOR/users/dcatlett/ampData_18sV9/Raw_data/",Lib10)


# designate a path to store all your plots of processing results. 
plotz_path1 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib1,"/processing_results/")
plotz_path2 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib2,"/processing_results/")
plotz_path3 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib3,"/processing_results/")
plotz_path4 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib4,"/processing_results/")
plotz_path5 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib5,"/processing_results/")
plotz_path6 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib6,"/processing_results/")
plotz_path7 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib7,"/processing_results/")
plotz_path8 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib8,"/processing_results/")
plotz_path9 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib9,"/processing_results/")
plotz_path10 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib10,"/processing_results/")


# create directories for storing plots of each library:
dir.create(plotz_path1, recursive = TRUE)
dir.create(plotz_path2, recursive = TRUE)
dir.create(plotz_path3, recursive = TRUE)
dir.create(plotz_path4, recursive = TRUE)
dir.create(plotz_path5, recursive = TRUE)
dir.create(plotz_path6, recursive = TRUE)
dir.create(plotz_path7, recursive = TRUE)
dir.create(plotz_path8, recursive = TRUE)
dir.create(plotz_path9, recursive = TRUE)
dir.create(plotz_path10, recursive = TRUE)


#### Workflow for Miseq data
### ----------------------------
## Step 1:
# Filter and Trim sequences with dada2_ML_filtering_Nov2018_dc.R

### paths for filtering script - change these:
## Where your raw reads are:
pathF1 <- path1 # see above^
pathR1 <- path1 # see above^
pathF2 <- path2 # see above^
pathR2 <- path2 # see above^
pathF3 <- path3 # see above^
pathR3 <- path3 # see above^
pathF4 <- path4 # see above^
pathR4 <- path4 # see above^
pathF5 <- path5 # see above^
pathR5 <- path5 # see above^
pathF6 <- path6 # see above^
pathR6 <- path6 # see above^
pathF7 <- path7 # see above^
pathR7 <- path7 # see above^
pathF8 <- path8 # see above^
pathR8 <- path8 # see above^
pathF9 <- path9 # see above^
pathR9 <- path9 # see above^
pathF10 <- path10 # see above^
pathR10 <- path10 # see above^

# Check Quality Scores first here:
source("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/dada2_ML_inspect_Qscores_Mar20.R")

## Where you want the filtered reads to go:
filtMethod <- "" # this can be empty if you want, but if not make it an indication of filter parameters you used.
filtpathF1 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib1,"filtered_reads",filtMethod) # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR1 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib1,"filtered_reads",filtMethod) # ...
filtpathF2 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib2,"filtered_reads",filtMethod) # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR2 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib2,"filtered_reads",filtMethod) # ...
filtpathF3 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib3,"filtered_reads",filtMethod) # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR3 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib3,"filtered_reads",filtMethod) # ...
filtpathF4 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib4,"filtered_reads",filtMethod) # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR4 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib4,"filtered_reads",filtMethod) # ...
filtpathF5 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib5,"filtered_reads",filtMethod) # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR5 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib5,"filtered_reads",filtMethod) # ...
filtpathF6 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib6,"filtered_reads",filtMethod) # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR6 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib6,"filtered_reads",filtMethod) # ...
filtpathF7 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib7,"filtered_reads",filtMethod) # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR7 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib7,"filtered_reads",filtMethod) # ...
filtpathF8 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib8,"filtered_reads",filtMethod) # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR8 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib8,"filtered_reads",filtMethod) # ...
filtpathF9 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib9,"filtered_reads",filtMethod) # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR9 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib9,"filtered_reads",filtMethod) # ...
filtpathF10 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib10,"filtered_reads",filtMethod) # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR10 <- file.path("/home/dcatlett/amp_data_processing/dada2_ML_Mar20",Lib10,"filtered_reads",filtMethod) # ...

## Where you want the readsin/out .csv file to go:
pathout1 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib1,"/processing_results/filtering/",filtMethod)
pathout2 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib2,"/processing_results/filtering/",filtMethod)
pathout3 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib3,"/processing_results/filtering/",filtMethod)
pathout4 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib4,"/processing_results/filtering/",filtMethod)
pathout5 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib5,"/processing_results/filtering/",filtMethod)
pathout6 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib6,"/processing_results/filtering/",filtMethod)
pathout7 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib7,"/processing_results/filtering/",filtMethod)
pathout8 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib8,"/processing_results/filtering/",filtMethod)
pathout9 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib9,"/processing_results/filtering/",filtMethod)
pathout10 <- paste0("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/",Lib10,"/processing_results/filtering/",filtMethod)

source("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/dada2_ML_filtering_Mar20_dc.R")

### ------------------------------
## Step 2:
## Infer Sequence Variants with dada2_big_data_infer_seq_variants_Jul2018_dc.R

### Paths for this script
## Designate the path you want to use to grab  the filtered reads. If you're not going straight through with
## the data set you filtered above, specify where you want to draw them from below. Also specify the filt method
# filtMethod <- ""
# filtpathF <- ""
# filtpathR <- ""

# same thing as filt method, just keeping track of parameter changes. This get's added onto the filtMethod
# in the directory name:
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

dir.create(pathout_ISV1)
dir.create(pathout_ISV2)
dir.create(pathout_ISV3)
dir.create(pathout_ISV4)
dir.create(pathout_ISV5)
dir.create(pathout_ISV6)
dir.create(pathout_ISV7)
dir.create(pathout_ISV8)
dir.create(pathout_ISV9)
dir.create(pathout_ISV10)

source("/home/dcatlett/amp_data_processing/dada2_ML_Mar20/dada2_ML_infer_seq_variants_Mar20_dc.R")

### If you made it this far with no breakage, continue to part 2 to combine libraries, remove chimeras, etc...
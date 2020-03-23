### dada2 pipeline 
## quality trimming and filtering raw MiSeq reads
## written by D. Catlett Mar 2020

## At present, the readsin/out file only includes read1. This is because right now we have it set up to filter R1
## and R2 simultaneously (still independently), and it throws out reads where both R1 and R2 don't pass the filter

## Library 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# File parsing
fastqFs <- sort(list.files(pathF1, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR1, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#
out <- filterAndTrim(fwd=file.path(pathF1, fastqFs), filt=file.path(filtpathF1, fastqFs),
                     rev=file.path(pathR1, fastqRs), filt.rev=file.path(filtpathR1, fastqRs),
                     truncLen=c(140,120), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

dir.create(pathout1)
write.csv(out, paste0(pathout1,"/readsin_readsout.csv"))

## Library 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# File parsing
fastqFs <- sort(list.files(pathF2, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR2, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#
out <- filterAndTrim(fwd=file.path(pathF2, fastqFs), filt=file.path(filtpathF2, fastqFs),
                     rev=file.path(pathR2, fastqRs), filt.rev=file.path(filtpathR2, fastqRs),
                     truncLen=c(140,120), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

dir.create(pathout2)
write.csv(out, paste0(pathout2,"/readsin_readsout.csv"))

## Library 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# File parsing
fastqFs <- sort(list.files(pathF3, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR3, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#
out <- filterAndTrim(fwd=file.path(pathF3, fastqFs), filt=file.path(filtpathF3, fastqFs),
                     rev=file.path(pathR3, fastqRs), filt.rev=file.path(filtpathR3, fastqRs),
                     truncLen=c(140,120), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

dir.create(pathout3)
write.csv(out, paste0(pathout3,"/readsin_readsout.csv"))

## Library 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# File parsing
fastqFs <- sort(list.files(pathF4, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR4, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#
out <- filterAndTrim(fwd=file.path(pathF4, fastqFs), filt=file.path(filtpathF4, fastqFs),
                     rev=file.path(pathR4, fastqRs), filt.rev=file.path(filtpathR4, fastqRs),
                     truncLen=c(140,120), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

dir.create(pathout4)
write.csv(out, paste0(pathout4,"/readsin_readsout.csv"))

## Library 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# File parsing
fastqFs <- sort(list.files(pathF5, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR5, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#
out <- filterAndTrim(fwd=file.path(pathF5, fastqFs), filt=file.path(filtpathF5, fastqFs),
                     rev=file.path(pathR5, fastqRs), filt.rev=file.path(filtpathR5, fastqRs),
                     truncLen=c(140,120), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

dir.create(pathout5)
write.csv(out, paste0(pathout5,"/readsin_readsout.csv"))

## Library 6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# File parsing
fastqFs <- sort(list.files(pathF6, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR6, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#
out <- filterAndTrim(fwd=file.path(pathF6, fastqFs), filt=file.path(filtpathF6, fastqFs),
                     rev=file.path(pathR6, fastqRs), filt.rev=file.path(filtpathR6, fastqRs),
                     truncLen=c(140,120), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

dir.create(pathout6)
write.csv(out, paste0(pathout6,"/readsin_readsout.csv"))

## Library 7 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# File parsing
fastqFs <- sort(list.files(pathF7, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR7, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#
out <- filterAndTrim(fwd=file.path(pathF7, fastqFs), filt=file.path(filtpathF7, fastqFs),
                     rev=file.path(pathR7, fastqRs), filt.rev=file.path(filtpathR7, fastqRs),
                     truncLen=c(140,120), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

dir.create(pathout7)
write.csv(out, paste0(pathout7,"/readsin_readsout.csv"))

## Library 8 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# File parsing
fastqFs <- sort(list.files(pathF8, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR8, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#
out <- filterAndTrim(fwd=file.path(pathF8, fastqFs), filt=file.path(filtpathF8, fastqFs),
                     rev=file.path(pathR8, fastqRs), filt.rev=file.path(filtpathR8, fastqRs),
                     truncLen=c(140,120), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

dir.create(pathout8)
write.csv(out, paste0(pathout8,"/readsin_readsout.csv"))

## Library 9 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# File parsing
fastqFs <- sort(list.files(pathF9, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR9, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#
out <- filterAndTrim(fwd=file.path(pathF9, fastqFs), filt=file.path(filtpathF9, fastqFs),
                     rev=file.path(pathR9, fastqRs), filt.rev=file.path(filtpathR9, fastqRs),
                     truncLen=c(140,120), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

dir.create(pathout9)
write.csv(out, paste0(pathout9,"/readsin_readsout.csv"))

## Library 10 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# File parsing
fastqFs <- sort(list.files(pathF10, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR10, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#
out <- filterAndTrim(fwd=file.path(pathF10, fastqFs), filt=file.path(filtpathF10, fastqFs),
                     rev=file.path(pathR10, fastqRs), filt.rev=file.path(filtpathR10, fastqRs),
                     truncLen=c(140,120), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

dir.create(pathout10)
write.csv(out, paste0(pathout10,"/readsin_readsout.csv"))


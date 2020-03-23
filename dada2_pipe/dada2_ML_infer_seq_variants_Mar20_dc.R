### This is part of Catlett's bioinformatic pipeline
## It does ASV inference, dereplication, and merging of paired reads for each library.

### Library 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# File parsing
filtFs <- list.files(filtpathF1, pattern="_R1_001.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR1, pattern="_R2_001.fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

## learn error rates:
errF <- learnErrors(filtFs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# # visualize and save error rates
# plotErrors(errF, nominalQ=TRUE)
plot.FerrMdl <- plotErrors(errF, nominalQ=TRUE)
ggsave(paste0(plotz_path1,"F_errMdl.pdf"), plot.FerrMdl, device="pdf")

# plotErrors(errR, nominalQ=TRUE)
plot.RerrMdl <- plotErrors(errR, nominalQ=TRUE)
ggsave(paste0(plotz_path1,"R_errMdl.pdf"), plot.RerrMdl, device="pdf")

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang = TRUE)
  mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

# keep track of reads merged and tabled:
Nreads.mergedAndTabled <- rowSums(seqtab)
write.csv(Nreads.mergedAndTabled, paste0(pathout_ISV1,"readsMergedAndTabled.csv"))

# Inspect distribution of sequence lengths
ampSizes <- table(nchar(getSequences(seqtab)))

saveRDS(seqtab, paste0(pathout_ISV1,"seqtab.rds")) # CHANGE ME to where you want sequence table saved

### Library 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# File parsing
filtFs <- list.files(filtpathF2, pattern="_R1_001.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR2, pattern="_R2_001.fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

## learn error rates:
errF <- learnErrors(filtFs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# # visualize and save error rates
plot.FerrMdl <- plotErrors(errF, nominalQ=TRUE)
ggsave(paste0(plotz_path2,"F_errMdl.pdf"), plot.FerrMdl, device="pdf")

plot.RerrMdl <- plotErrors(errR, nominalQ=TRUE)
ggsave(paste0(plotz_path2,"R_errMdl.pdf"), plot.RerrMdl, device="pdf")

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang = TRUE)
  mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

# keep track of reads merged and tabled:
Nreads.mergedAndTabled <- rowSums(seqtab)
write.csv(Nreads.mergedAndTabled, paste0(pathout_ISV2,"readsMergedAndTabled.csv"))

# Inspect distribution of sequence lengths
ampSizes <- table(nchar(getSequences(seqtab)))

saveRDS(seqtab, paste0(pathout_ISV2,"seqtab.rds")) # CHANGE ME to where you want sequence table saved


### Library 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# File parsing
filtFs <- list.files(filtpathF3, pattern="_R1_001.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR3, pattern="_R2_001.fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

## learn error rates:
errF <- learnErrors(filtFs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# # visualize and save error rates
# plotErrors(errF, nominalQ=TRUE)
plot.FerrMdl <- plotErrors(errF, nominalQ=TRUE)
ggsave(paste0(plotz_path3,"F_errMdl.pdf"), plot.FerrMdl, device="pdf")

# plotErrors(errR, nominalQ=TRUE)
plot.RerrMdl <- plotErrors(errR, nominalQ=TRUE)
ggsave(paste0(plotz_path3,"R_errMdl.pdf"), plot.RerrMdl, device="pdf")

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang = TRUE)
  mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

# keep track of reads merged and tabled:
Nreads.mergedAndTabled <- rowSums(seqtab)
write.csv(Nreads.mergedAndTabled, paste0(pathout_ISV3,"readsMergedAndTabled.csv"))

# Inspect distribution of sequence lengths
ampSizes <- table(nchar(getSequences(seqtab)))

saveRDS(seqtab, paste0(pathout_ISV3,"seqtab.rds")) # CHANGE ME to where you want sequence table saved

### Library 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# File parsing
filtFs <- list.files(filtpathF4, pattern="_R1_001.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR4, pattern="_R2_001.fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

## learn error rates:
errF <- learnErrors(filtFs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# # visualize and save error rates
# plotErrors(errF, nominalQ=TRUE)
plot.FerrMdl <- plotErrors(errF, nominalQ=TRUE)
ggsave(paste0(plotz_path4,"F_errMdl.pdf"), plot.FerrMdl, device="pdf")

# plotErrors(errR, nominalQ=TRUE)
plot.RerrMdl <- plotErrors(errR, nominalQ=TRUE)
ggsave(paste0(plotz_path4,"R_errMdl.pdf"), plot.RerrMdl, device="pdf")

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang = TRUE)
  mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

# keep track of reads merged and tabled:
Nreads.mergedAndTabled <- rowSums(seqtab)
write.csv(Nreads.mergedAndTabled, paste0(pathout_ISV4,"readsMergedAndTabled.csv"))

# Inspect distribution of sequence lengths
ampSizes <- table(nchar(getSequences(seqtab)))

saveRDS(seqtab, paste0(pathout_ISV4,"seqtab.rds")) # CHANGE ME to where you want sequence table saved

### Library 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# File parsing
filtFs <- list.files(filtpathF5, pattern="_R1_001.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR5, pattern="_R2_001.fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

## learn error rates:
errF <- learnErrors(filtFs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# # visualize and save error rates
# plotErrors(errF, nominalQ=TRUE)
plot.FerrMdl <- plotErrors(errF, nominalQ=TRUE)
ggsave(paste0(plotz_path5,"F_errMdl.pdf"), plot.FerrMdl, device="pdf")

# plotErrors(errR, nominalQ=TRUE)
plot.RerrMdl <- plotErrors(errR, nominalQ=TRUE)
ggsave(paste0(plotz_path5,"R_errMdl.pdf"), plot.RerrMdl, device="pdf")

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang = TRUE)
  mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

# keep track of reads merged and tabled:
Nreads.mergedAndTabled <- rowSums(seqtab)
write.csv(Nreads.mergedAndTabled, paste0(pathout_ISV5,"readsMergedAndTabled.csv"))

# Inspect distribution of sequence lengths
ampSizes <- table(nchar(getSequences(seqtab)))

saveRDS(seqtab, paste0(pathout_ISV5,"seqtab.rds")) # CHANGE ME to where you want sequence table saved

### Library 6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# File parsing
filtFs <- list.files(filtpathF6, pattern="_R1_001.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR6, pattern="_R2_001.fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

## learn error rates:
errF <- learnErrors(filtFs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# # visualize and save error rates
# plotErrors(errF, nominalQ=TRUE)
plot.FerrMdl <- plotErrors(errF, nominalQ=TRUE)
ggsave(paste0(plotz_path6,"F_errMdl.pdf"), plot.FerrMdl, device="pdf")

# plotErrors(errR, nominalQ=TRUE)
plot.RerrMdl <- plotErrors(errR, nominalQ=TRUE)
ggsave(paste0(plotz_path6,"R_errMdl.pdf"), plot.RerrMdl, device="pdf")

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang = TRUE)
  mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

# keep track of reads merged and tabled:
Nreads.mergedAndTabled <- rowSums(seqtab)
write.csv(Nreads.mergedAndTabled, paste0(pathout_ISV6,"readsMergedAndTabled.csv"))

# Inspect distribution of sequence lengths
ampSizes <- table(nchar(getSequences(seqtab)))

saveRDS(seqtab, paste0(pathout_ISV6,"seqtab.rds")) # CHANGE ME to where you want sequence table saved

### Library 7 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# File parsing
filtFs <- list.files(filtpathF7, pattern="_R1_001.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR7, pattern="_R2_001.fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

## learn error rates:
errF <- learnErrors(filtFs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# # visualize and save error rates
# plotErrors(errF, nominalQ=TRUE)
plot.FerrMdl <- plotErrors(errF, nominalQ=TRUE)
ggsave(paste0(plotz_path7,"F_errMdl.pdf"), plot.FerrMdl, device="pdf")

# plotErrors(errR, nominalQ=TRUE)
plot.RerrMdl <- plotErrors(errR, nominalQ=TRUE)
ggsave(paste0(plotz_path7,"R_errMdl.pdf"), plot.RerrMdl, device="pdf")

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang = TRUE)
  mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

# keep track of reads merged and tabled:
Nreads.mergedAndTabled <- rowSums(seqtab)
write.csv(Nreads.mergedAndTabled, paste0(pathout_ISV7,"readsMergedAndTabled.csv"))

# Inspect distribution of sequence lengths
ampSizes <- table(nchar(getSequences(seqtab)))

saveRDS(seqtab, paste0(pathout_ISV7,"seqtab.rds")) # CHANGE ME to where you want sequence table saved

### Library 8 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# File parsing
filtFs <- list.files(filtpathF8, pattern="_R1_001.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR8, pattern="_R2_001.fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

## learn error rates:
errF <- learnErrors(filtFs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# # visualize and save error rates
# plotErrors(errF, nominalQ=TRUE)
plot.FerrMdl <- plotErrors(errF, nominalQ=TRUE)
ggsave(paste0(plotz_path8,"F_errMdl.pdf"), plot.FerrMdl, device="pdf")

# plotErrors(errR, nominalQ=TRUE)
plot.RerrMdl <- plotErrors(errR, nominalQ=TRUE)
ggsave(paste0(plotz_path8,"R_errMdl.pdf"), plot.RerrMdl, device="pdf")

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang = TRUE)
  mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

# keep track of reads merged and tabled:
Nreads.mergedAndTabled <- rowSums(seqtab)
write.csv(Nreads.mergedAndTabled, paste0(pathout_ISV8,"readsMergedAndTabled.csv"))

# Inspect distribution of sequence lengths
ampSizes <- table(nchar(getSequences(seqtab)))

saveRDS(seqtab, paste0(pathout_ISV8,"seqtab.rds")) # CHANGE ME to where you want sequence table saved

### Library 9 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# File parsing
filtFs <- list.files(filtpathF9, pattern="_R1_001.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR9, pattern="_R2_001.fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

## learn error rates:
errF <- learnErrors(filtFs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# # visualize and save error rates
# plotErrors(errF, nominalQ=TRUE)
plot.FerrMdl <- plotErrors(errF, nominalQ=TRUE)
ggsave(paste0(plotz_path9,"F_errMdl.pdf"), plot.FerrMdl, device="pdf")

# plotErrors(errR, nominalQ=TRUE)
plot.RerrMdl <- plotErrors(errR, nominalQ=TRUE)
ggsave(paste0(plotz_path9,"R_errMdl.pdf"), plot.RerrMdl, device="pdf")

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
    cat("Processing:", sam, "\n")
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang = TRUE)
    mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

# keep track of reads merged and tabled:
Nreads.mergedAndTabled <- rowSums(seqtab)
write.csv(Nreads.mergedAndTabled, paste0(pathout_ISV9,"readsMergedAndTabled.csv"))

# Inspect distribution of sequence lengths
ampSizes <- table(nchar(getSequences(seqtab)))

saveRDS(seqtab, paste0(pathout_ISV9,"seqtab.rds")) # CHANGE ME to where you want sequence table saved


### Library 10 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# File parsing
filtFs <- list.files(filtpathF10, pattern="_R1_001.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR10, pattern="_R2_001.fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

## learn error rates:
errF <- learnErrors(filtFs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# # visualize and save error rates
# plotErrors(errF, nominalQ=TRUE)
plot.FerrMdl <- plotErrors(errF, nominalQ=TRUE)
ggsave(paste0(plotz_path10,"F_errMdl.pdf"), plot.FerrMdl, device="pdf")

# plotErrors(errR, nominalQ=TRUE)
plot.RerrMdl <- plotErrors(errR, nominalQ=TRUE)
ggsave(paste0(plotz_path10,"R_errMdl.pdf"), plot.RerrMdl, device="pdf")

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
    cat("Processing:", sam, "\n")
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang = TRUE)
    mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

# keep track of reads merged and tabled:
Nreads.mergedAndTabled <- rowSums(seqtab)
write.csv(Nreads.mergedAndTabled, paste0(pathout_ISV10,"readsMergedAndTabled.csv"))

# Inspect distribution of sequence lengths
ampSizes <- table(nchar(getSequences(seqtab)))

saveRDS(seqtab, paste0(pathout_ISV10,"seqtab.rds")) # CHANGE ME to where you want sequence table saved



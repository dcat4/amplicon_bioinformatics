### This script written by D. Catlett, March 2020. 
### Just plots a quick check of a random sub-sample of Q-scores from F and R reads, as outlined in dada2 tutorial. 

# Path's are set up in Dylan's pipeline progression.

#### Library 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## grab fastq's
fastqFs <- sort(list.files(pathF1, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR1, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# plot Quality scores; may need to subset if you're library is heavily multiplexed
# saves to a pdf file for review
plot.Rquals <- plotQualityProfile(paste0(pathR1,"/",fastqRs))
ggsave(paste0(plotz_path1,"Rqualplot.pdf"), plot.Rquals, device="pdf")

plot.Fquals <- plotQualityProfile(paste0(pathF1,"/",fastqFs))
ggsave(paste0(plotz_path1,"Fqualplot.pdf"), plot.Fquals, device="pdf")


#### Library 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## grab fastq's
fastqFs <- sort(list.files(pathF2, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR2, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# plot Quality scores; may need to subset if you're library is heavily multiplexed
# saves to a pdf file for review
plot.Rquals <- plotQualityProfile(paste0(pathR2,"/",fastqRs))
ggsave(paste0(plotz_path2,"Rqualplot.pdf"), plot.Rquals, device="pdf")

plot.Fquals <- plotQualityProfile(paste0(pathF2,"/",fastqFs))
ggsave(paste0(plotz_path2,"Fqualplot.pdf"), plot.Fquals, device="pdf")

#### Library 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## grab fastq's
fastqFs <- sort(list.files(pathF3, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR3, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# plot Quality scores; may need to subset if you're library is heavily multiplexed
# saves to a pdf file for review
plot.Rquals <- plotQualityProfile(paste0(pathR3,"/",fastqRs))
ggsave(paste0(plotz_path3,"Rqualplot.pdf"), plot.Rquals, device="pdf")

plot.Fquals <- plotQualityProfile(paste0(pathF3,"/",fastqFs))
ggsave(paste0(plotz_path3,"Fqualplot.pdf"), plot.Fquals, device="pdf")

#### Library 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## grab fastq's
fastqFs <- sort(list.files(pathF4, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR4, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# plot Quality scores; may need to subset if you're library is heavily multiplexed
# saves to a pdf file for review
plot.Rquals <- plotQualityProfile(paste0(pathR4,"/",fastqRs))
ggsave(paste0(plotz_path4,"Rqualplot.pdf"), plot.Rquals, device="pdf")

plot.Fquals <- plotQualityProfile(paste0(pathF4,"/",fastqFs))
ggsave(paste0(plotz_path4,"Fqualplot.pdf"), plot.Fquals, device="pdf")

#### Library 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## grab fastq's
fastqFs <- sort(list.files(pathF5, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR5, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# plot Quality scores; may need to subset if you're library is heavily multiplexed
# saves to a pdf file for review
plot.Rquals <- plotQualityProfile(paste0(pathR5,"/",fastqRs))
ggsave(paste0(plotz_path5,"Rqualplot.pdf"), plot.Rquals, device="pdf")

plot.Fquals <- plotQualityProfile(paste0(pathF5,"/",fastqFs))
ggsave(paste0(plotz_path5,"Fqualplot.pdf"), plot.Fquals, device="pdf")

#### Library 6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## grab fastq's
fastqFs <- sort(list.files(pathF6, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR6, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# plot Quality scores; may need to subset if you're library is heavily multiplexed
# saves to a pdf file for review
plot.Rquals <- plotQualityProfile(paste0(pathR6,"/",fastqRs))
ggsave(paste0(plotz_path6,"Rqualplot.pdf"), plot.Rquals, device="pdf")

plot.Fquals <- plotQualityProfile(paste0(pathF6,"/",fastqFs))
ggsave(paste0(plotz_path6,"Fqualplot.pdf"), plot.Fquals, device="pdf")

#### Library 7 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## grab fastq's
fastqFs <- sort(list.files(pathF7, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR7, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# plot Quality scores; may need to subset if you're library is heavily multiplexed
# saves to a pdf file for review
plot.Rquals <- plotQualityProfile(paste0(pathR7,"/",fastqRs))
ggsave(paste0(plotz_path7,"Rqualplot.pdf"), plot.Rquals, device="pdf")

plot.Fquals <- plotQualityProfile(paste0(pathF7,"/",fastqFs))
ggsave(paste0(plotz_path7,"Fqualplot.pdf"), plot.Fquals, device="pdf")

#### Library 8 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## grab fastq's
fastqFs <- sort(list.files(pathF8, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR8, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# plot Quality scores; may need to subset if you're library is heavily multiplexed
# saves to a pdf file for review
plot.Rquals <- plotQualityProfile(paste0(pathR8,"/",fastqRs))
ggsave(paste0(plotz_path8,"Rqualplot.pdf"), plot.Rquals, device="pdf")

plot.Fquals <- plotQualityProfile(paste0(pathF8,"/",fastqFs))
ggsave(paste0(plotz_path8,"Fqualplot.pdf"), plot.Fquals, device="pdf")

#### Library 9 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## grab fastq's
fastqFs <- sort(list.files(pathF9, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR9, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# plot Quality scores; may need to subset if you're library is heavily multiplexed
# saves to a pdf file for review
plot.Rquals <- plotQualityProfile(paste0(pathR9,"/",fastqRs))
ggsave(paste0(plotz_path9,"Rqualplot.pdf"), plot.Rquals, device="pdf")

plot.Fquals <- plotQualityProfile(paste0(pathF9,"/",fastqFs))
ggsave(paste0(plotz_path9,"Fqualplot.pdf"), plot.Fquals, device="pdf")

#### Library 10 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## grab fastq's
fastqFs <- sort(list.files(pathF10, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR10, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# plot Quality scores; may need to subset if you're library is heavily multiplexed
# saves to a pdf file for review
plot.Rquals <- plotQualityProfile(paste0(pathR10,"/",fastqRs))
ggsave(paste0(plotz_path10,"Rqualplot.pdf"), plot.Rquals, device="pdf")

plot.Fquals <- plotQualityProfile(paste0(pathF10,"/",fastqFs))
ggsave(paste0(plotz_path10,"Fqualplot.pdf"), plot.Fquals, device="pdf")



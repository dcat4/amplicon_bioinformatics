# this script mines pr2 and silva databases to assemble a large test data set of V4 sequences...

# this assembles v4 amplicons from pr2 and silva independtently 
# does not do expected taxnomies b/c R is too slow -- 
# instead gonna save csv files of v4 asvs and try in matlab

rm(list=ls())
setwd("~/Documents/R/amplicon_bioinformatics/mock_analysis/")

library("DECIPHER")
library("Biostrings")
library("stringr")

fastaFile <- readDNAStringSet("~/Documents/R/pr2_version_4.12.0_18S_dada2.fasta")
seq_name = names(fastaFile)
# assign local accession codes to each sequence in silva for later by appending to seq_name
seq_name2 <- str_replace_all(seq_name, ";", "|")
for (i in 1:length(seq_name)) {
  seq_name[i] <- paste0("pr2.",toString(i),"|",seq_name2[i])
}
sequence = paste(fastaFile)
pr2 <- data.frame(seq_name, sequence, stringsAsFactors = FALSE)

fastaFile <- readDNAStringSet("~/Documents/R/silva_nr_v138_train_set.fa")
seq_name = names(fastaFile)
# assign local accession codes to each sequence in silva for later by appending to seq_name
seq_name2 <- str_replace_all(seq_name, ";", "|")
for (i in 1:length(seq_name)) {
  seq_name[i] <- paste0("silva",toString(i),"|",seq_name2[i])
}
sequence = paste(fastaFile)
silva <- data.frame(seq_name, sequence, stringsAsFactors = FALSE)
silvaeuk <- silva[which(str_detect(silva$seq_name,"Eukaryota")),]
hist(str_length(silvaeuk$sequence))

silva <- silvaeuk
# clear out unnecessary stuff:
rm("fastaFile", "seq_name", "seq_name2", "sequence")

# remove sequences in silva that are found in pr2 to make sure your ref data is all unique (this favors pr2)
ii <- which(silva$sequence %in% pr2$sequence)
i2 <- which(pr2$sequence %in% silva$sequence)

# isolate v4, and v4-5 regions
# adapted from here: https://github.com/pr2database/pr2-primers/blob/master/PR2%20Primers%20pr2_match.R
pp <- read.csv("primer_sets.csv", stringsAsFactors = FALSE)
ur <- unique(pp$X18s_region)

pr2.2 <- pr2$sequence
names(pr2.2) <- pr2$seq_name
pr2 <- DNAStringSet(pr2.2)

silva.2 <- silva$sequence
names(silva.2) <- silva$seq_name
silva <- DNAStringSet(silva.2)

mm <- 2 # max mismatches for finding primer hits

pr2.prime <- vector(mode = "list", length = length(ur)*2)
silva.prime <- vector(mode = "list", length = length(ur)*2)
counter <- 1
for (i in 1:length(ur)) {
  fwd <-  DNAString(pp$sequence[min(which(pp$X18s_region == ur[i]))])
  rev <-  DNAString(pp$sequence[max(which(pp$X18s_region == ur[i]))])
  rev <-  reverseComplement(rev)
  
  fwd.pos.pr2 <- vmatchPattern(fwd, pr2, max.mismatch=mm, min.mismatch=0, with.indels=FALSE, fixed=FALSE, algorithm="auto")
  fwd.pos.silva <- vmatchPattern(fwd, silva, max.mismatch=mm, min.mismatch=0, with.indels=FALSE, fixed=FALSE, algorithm="auto")
  rev.pos.pr2 <- vmatchPattern(rev, pr2, max.mismatch=mm, min.mismatch=0, with.indels=FALSE, fixed=FALSE, algorithm="auto")
  rev.pos.silva <- vmatchPattern(rev, silva, max.mismatch=mm, min.mismatch=0, with.indels=FALSE, fixed=FALSE, algorithm="auto")
  
  # unlist the MIndex position
  fwd.pos.pr2 <- unlist(fwd.pos.pr2)
  fwd.pos.silva <- unlist(fwd.pos.silva)
  rev.pos.pr2 <- unlist(rev.pos.pr2)
  rev.pos.silva <- unlist(rev.pos.silva)
  
  # Create data frames to store primer match data
  # pr2:
  fhits.pr2 <- data.frame(namer = fwd.pos.pr2@NAMES, 
                          fpos.strt = fwd.pos.pr2@start, 
                          fpos.end = fwd.pos.pr2@start + (length(fwd)-1),
                          stringsAsFactors = FALSE)
  rhits.pr2 <- data.frame(namer = rev.pos.pr2@NAMES, 
                          rpos.strt = rev.pos.pr2@start,
                          rpos.end = rev.pos.pr2@start + (length(rev)-1),
                          stringsAsFactors = FALSE)
  # silva:
  fhits.silva <- data.frame(namer = fwd.pos.silva@NAMES, 
                            fpos.strt = fwd.pos.silva@start, 
                            fpos.end = fwd.pos.silva@start + (length(fwd)-1),
                            stringsAsFactors = FALSE)
  rhits.silva <- data.frame(namer = rev.pos.silva@NAMES, 
                            rpos.strt = rev.pos.silva@start,
                            rpos.end = rev.pos.silva@start + (length(rev)-1),
                            stringsAsFactors = FALSE)
  # store primer hit data in a list:
  pr2.prime[[counter]] <- fhits.pr2
  names(pr2.prime)[counter] <- paste0("fwd.",ur[i])
  silva.prime[[counter]] <- fhits.silva
  names(silva.prime)[counter] <- paste0("fwd.",ur[i])
  counter <- counter + 1
  pr2.prime[[counter]] <- rhits.pr2
  names(pr2.prime)[counter] <- paste0("rev.",ur[i])
  silva.prime[[counter]] <- rhits.silva
  names(silva.prime)[counter] <- paste0("rev.",ur[i])
  counter <- counter + 1
}

# use primer hit results to populate dataframes with accession numbers/taxonomies/ASVs that can be used in a mock data set
hv <- "v4"
fsub <- pr2.prime[[which(str_detect(names(pr2.prime), paste0("fwd.", hv)) & !str_detect(names(pr2.prime), "5"))]]
fsub <- fsub[which(!is.na(fsub$namer)),]
rsub <- pr2.prime[[which(str_detect(names(pr2.prime), paste0("rev.", hv)) & !str_detect(names(pr2.prime), "5"))]]
rsub <- rsub[which(!is.na(rsub$namer)),]
pr2v <- as.character(pr2) # a vector of the reference sequnences...

bb <- intersect(fsub$namer, rsub$namer)
fsub <- fsub[fsub$namer %in% bb,]
rsub <- rsub[rsub$namer %in% bb,]
un <- unique(fsub$namer)
nav <- rep(NA, times = length(un))
v4.pr2 <- data.frame(acc.tax = un, db = rep("pr2", times = length(un)), start = nav, stop = nav, asv.len = nav, asv = nav, stringsAsFactors = FALSE)
bloop <- c()
for (i in 1:length(un)) {
  fhits <- fsub[fsub$namer == un[i],]
  rhits <- rsub[rsub$namer == un[i],]
  strt <- fhits$fpos.end + 1
  fin <- rhits$rpos.strt - 1
  
  if (length(strt) == 1 && length(fin) == 1) {
    
    if (fin <= strt) {
      bloop <- append(bloop, i)
      v4.pr2[i, c("start","stop","asv.len")] <- c(strt, fin, length(strt:fin))
    } else {
      v4.pr2[i, c("start","stop","asv.len")] <- c(strt, fin, length(strt:fin))
      ss <- pr2v[un[i]]
      v4.pr2[i, "asv"] <- str_sub(ss, start = strt, end = fin)
    }
    
  } else if (length(strt) == 1 && length(fin) > 1) {
    # account for multiple reverse hits
    # use the amplicon who's size is closer to target of 380 --> see Bradley and pr2 github page
    x <- fin-strt
    eh <- which(x > 0 & abs(x - 380) == min(abs(x - 380)))
    fin <- fin[eh]
    v4.pr2[i, c("start","stop","asv.len")] <- c(strt, fin, length(strt:fin))
    ss <- pr2v[un[i]]
    v4.pr2[i, "asv"] <- str_sub(ss, start = strt, end = fin)
    
  } else if (length(strt) > 1 && length(fin) == 1) {
    x <- fin-strt
    eh <- which(x > 0 & abs(x - 380) == min(abs(x - 380)))
    strt <- strt[eh]
    v4.pr2[i, c("start","stop","asv.len")] <- c(strt, fin, length(strt:fin))
    ss <- pr2v[un[i]]
    v4.pr2[i, "asv"] <- str_sub(ss, start = strt, end = fin)
    
  } else if (length(strt) > 1 && length(fin) > 1) {
    # account for multiple F + R hits
    eh <- sapply(strt, function(x) fin - x)
    ehF <- apply(eh, MARGIN = 2, function(x) which(x > 0 & abs(x - 380) == min(abs(x - 380))))
    ehR <- apply(eh, MARGIN = 1, function(x) which(x > 0 & abs(x - 380) == min(abs(x - 380))))
    if (length(ehF) != length(ehR)){
      error("issue in amplicon finding where both primers hit more than once")
    }
    # ehF and ehR are vectors that are = lengths
    # corresponding entries in the 2 vectors are pairs of strt/stop positions of target ASVs
    # if lenght 1 you can just use the single target ASV. if >length 1 just keep the the first one 
    if (length(ehF) == 1) {
      strt <- strt[ehF]
      fin <- fin[ehR]
      v4.pr2[i, c("start","stop","asv.len")] <- c(strt, fin, length(strt:fin))
      ss <- pr2v[un[i]]
      v4.pr2[i, "asv"] <- str_sub(ss, start = strt, end = fin)
    } else if (length(ehF) > 1) {
      strt <- strt[ehF[1]]
      fin <- fin[ehR[1]]
      v4.pr2[i, c("start","stop","asv.len")] <- c(strt, fin, length(strt:fin))
      ss <- pr2v[un[i]]
      v4.pr2[i, "asv"] <- str_sub(ss, start = strt, end = fin)
    }
  }
  
}
# remove bloop (where reverse hit further 5' than fwd) --> these are NA's so just rm NA's
v4.pr2 <- v4.pr2[which(!is.na(v4.pr2$asv)) , ]

### repeat above for silva-v4:
fsub <- silva.prime[[which(str_detect(names(silva.prime), paste0("fwd.", hv)) & !str_detect(names(silva.prime), "5"))]]
fsub <- fsub[which(!is.na(fsub$namer)),]
rsub <- silva.prime[[which(str_detect(names(silva.prime), paste0("rev.", hv)) & !str_detect(names(silva.prime), "5"))]]
rsub <- rsub[which(!is.na(rsub$namer)),]
silvav <- as.character(silva) # a vector of the reference sequnences...

bb <- intersect(fsub$namer, rsub$namer)
fsub <- fsub[fsub$namer %in% bb,]
rsub <- rsub[rsub$namer %in% bb,]
un <- unique(fsub$namer)
nav <- rep(NA, times = length(un))
v4.silva <- data.frame(acc.tax = un, db = rep("silva", times = length(un)), start = nav, stop = nav, asv.len = nav, asv = nav, stringsAsFactors = FALSE)
bloop <- c()
for (i in 1:length(un)) {
  fhits <- fsub[fsub$namer == un[i],]
  rhits <- rsub[rsub$namer == un[i],]
  strt <- fhits$fpos.end + 1
  fin <- rhits$rpos.strt - 1
  
  if (length(strt) == 1 && length(fin) == 1) {
    
    if (fin <= strt) {
      bloop <- append(bloop, i)
      v4.silva[i, c("start","stop","asv.len")] <- c(strt, fin, length(strt:fin))
    } else {
      v4.silva[i, c("start","stop","asv.len")] <- c(strt, fin, length(strt:fin))
      ss <- silvav[un[i]]
      v4.silva[i, "asv"] <- str_sub(ss, start = strt, end = fin)
    }
    
  } else if (length(strt) == 1 && length(fin) > 1) {
    # account for multiple reverse hits
    # use the amplicon who's size is closer to target of 380 --> see Bradley and pr2 github page
    x <- fin-strt
    eh <- which(x > 0 & abs(x - 380) == min(abs(x - 380)))
    fin <- fin[eh]
    v4.silva[i, c("start","stop","asv.len")] <- c(strt, fin, length(strt:fin))
    ss <- silvav[un[i]]
    v4.silva[i, "asv"] <- str_sub(ss, start = strt, end = fin)
    
  } else if (length(strt) > 1 && length(fin) == 1) {
    x <- fin-strt
    eh <- which(x > 0 & abs(x - 380) == min(abs(x - 380)))
    strt <- strt[eh]
    v4.silva[i, c("start","stop","asv.len")] <- c(strt, fin, length(strt:fin))
    ss <- silvav[un[i]]
    v4.silva[i, "asv"] <- str_sub(ss, start = strt, end = fin)
    
  } else if (length(strt) > 1 && length(fin) > 1) {
    # account for multiple F + R hits
    eh <- sapply(strt, function(x) fin - x)
    ehF <- apply(eh, MARGIN = 2, function(x) which(x > 0 & abs(x - 380) == min(abs(x - 380))))
    ehR <- apply(eh, MARGIN = 1, function(x) which(x > 0 & abs(x - 380) == min(abs(x - 380))))
    if (length(ehF) != length(ehR)){
      error("issue in amplicon finding where both primers hit more than once")
    }
    # ehF and ehR are vectors that are = lengths
    # corresponding entries in the 2 vectors are pairs of strt/stop positions of target ASVs
    # if lenght 1 you can just use the single target ASV. if >length 1 just keep the the first one 
    if (length(ehF) == 1) {
      strt <- strt[ehF]
      fin <- fin[ehR]
      v4.silva[i, c("start","stop","asv.len")] <- c(strt, fin, length(strt:fin))
      ss <- silvav[un[i]]
      v4.silva[i, "asv"] <- str_sub(ss, start = strt, end = fin)
    } else if (length(ehF) > 1) {
      strt <- strt[ehF[1]]
      fin <- fin[ehR[1]]
      v4.silva[i, c("start","stop","asv.len")] <- c(strt, fin, length(strt:fin))
      ss <- silvav[un[i]]
      v4.silva[i, "asv"] <- str_sub(ss, start = strt, end = fin)
    }
  }
  
}

# remove bloop (where reverse hit further 5' than fwd) --> these are NA's so just rm NA's
v4.silva <- v4.silva[which(!is.na(v4.silva$asv)) , ]

# got all the v4 shit done and stored in 4 datasets
# cleaning should probably happen -- deal w/ 
# 4. merging the 2 db's
# 5. check that hypervariable ASVs are still unique and trim to a common taxonomy if not (like LCA style)...
foooooooooook
# check for NA's in your 4 data sets:
if (any(is.na(v4.pr2)) || any(is.na(v4.silva))) {
  error("NA's present in one of your v4 datasets")
}

# look at ASV lengths from the 4 data sets:
boxplot(c(v4.pr2$asv.len, v4.silva$asv.len))
v4.pr2 <- v4.pr2[v4.pr2$asv.len > 90 & v4.pr2$asv.len < 180 ,]
v4.silva <- v4.silva[v4.silva$asv.len > 90 & v4.pr2$silva.len < 180 ,]
# check that it worked:
boxplot(c(v4.pr2$asv.len, v4.silva$asv.len))

# remove asv's that contain a degeneracy:
v4.pr2 <- v4.pr2[str_detect(v4.pr2$asv,"N",negate=TRUE) ,]
v4.silva <- v4.silva[str_detect(v4.silva$asv,"N",negate=TRUE) ,]

# check for multiple primer hits on the same ref:
if (length(which((duplicated(v4.pr2$acc.tax) | duplicated(v4.pr2$acc.tax, fromLast = TRUE)))) != 0) {
  error("need to deal w/ multiple hits in v4.pr2")
}

if (length(which((duplicated(v4.silva$acc.tax) | duplicated(v4.silva$acc.tax, fromLast = TRUE)))) != 0) {
  error("need to deal w/ multiple hits in v4.silva")
}

# write out your v4 ASVs to csv files so you can put the expected taxonomies together in matlab
# (otherwise it's mega-slow b/c R sux)
write.csv(v4.pr2, file = "mock_data/pr2_v4amps.csv")
write.csv(v4.silva, file = "mock_data/silva_v4amps.csv")



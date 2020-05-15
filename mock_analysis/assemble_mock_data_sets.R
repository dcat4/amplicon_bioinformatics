# this script mines pr2 and silva databases to assemble a large test data set of V9 and V4 sequences...

# figured out the issues w/ silva to some extent - reading in the fasta is removing T's (or U's). LOL!

# list of to-dos below for a final clean-up
# also need to make sure NA-containing rows are being removed from all the final data.frames

rm(list=ls())
setwd("~/Documents/R/amplicon_bioinformatics/mock_analysis/")
source("~/Documents/R/amplicon_bioinformatics/package_deal/all_of_it.R")

library("pr2database")
library("DECIPHER")
library("Biostrings")
library("stringr")

data("pr2")

# the fields are described here:
# https://pr2-database.org/documentation/pr2-fields/

# extract reference sequences (according to pr2 creator)
refz <- pr2[!is.na(pr2$reference_sequence),]
# this chex out:
# > unique(refz$reference_sequence)
# [1] 1
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
rm("pr2", "fastaFile", "seq_name", "sequence")

# pull sequences from refz, fix up names to match silva, and combine silva/pr2
# pr2 names should go 
refz <- subset(refz, select = c("pr2_accession","species","kingdom","supergroup","division","class","order","family","genus","sequence"))

df <- as.character(NA, ncol = 1, nrow = nrow(refz))
for (i in 1:nrow(refz)) {
  # xx <- paste0(refz$pr2_accession[i], " ", refz$kingdom[i], "; ", refz$supergroup[i], "; ", refz$division[i], "; ", refz$class[i], "; "
  #              , refz$order[i], "; ", refz$family[i], "; ", refz$genus[i], "; ", refz$species[i])
  xx <- paste0(refz$pr2_accession[i],"|", refz$kingdom[i],"|", refz$supergroup[i],"|", refz$division[i],"|", refz$class[i],"|",
               refz$order[i],"|", refz$family[i],"|", refz$genus[i], "|",refz$species[i])
  df[i] <- xx
}
# worked, hooray! now df is a character vector w/ pr2 accession and taxonomy. 

pr2 <- data.frame(df, stringsAsFactors = FALSE)
pr2$sequence <- refz$sequence
colnames(pr2) <- colnames(silva)

# remove sequences in silva that are found in pr2 to make sure your ref data is all unique (this favors pr2)
ii <- which(silva$sequence %in% pr2$sequence)
i2 <- which(pr2$sequence %in% silva$sequence)

# shockingly, there are no overlapping sequences between silva + pr2. there are a couple hundred duplicate sequences in each though.
# remove duplicate sequences from each:
pr2 <- pr2[!(duplicated(pr2$sequence) | duplicated(pr2$sequence, fromLast = TRUE)), ]
silva <- silva[!(duplicated(silva$sequence) | duplicated(silva$sequence, fromLast = TRUE)), ]

# isolate v9, v4, and v4-5 regions
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
hv <- "v9"
fsub <- pr2.prime[[which(str_detect(names(pr2.prime), paste0("fwd.", hv)))]]
fsub <- fsub[which(!is.na(fsub$namer)),]
rsub <- pr2.prime[[which(str_detect(names(pr2.prime), paste0("rev.", hv)))]]
rsub <- rsub[which(!is.na(rsub$namer)),]
pr2v <- as.character(pr2) # a vector of the reference sequnences...

bb <- setdiff(fsub$namer, rsub$namer)
bloop <- fsub[fsub$namer %in% bb,]
strt <- bloop$fpos.end + 1
ss <- pr2v[bloop$namer]
ll <- str_length(ss) - (strt+1)

v9.pr2.noR <- data.frame(acc.tax = bloop$namer, db = rep("pr2", times = length(bloop$namer)), start = strt, stop = str_length(ss), 
                         asv.len = ll, asv = rep(NA, times = length(bloop$namer)), stringsAsFactors = FALSE)
for (i in 1:nrow(v9.pr2.noR)) {
  v9.pr2.noR[i, "asv"] <- str_sub(ss[i], start = v9.pr2.noR[i, "start"], end = length(ss))
}

bb <- intersect(fsub$namer, rsub$namer)
fsub <- fsub[fsub$namer %in% bb,]
rsub <- rsub[rsub$namer %in% bb,]
un <- unique(fsub$namer)
nav <- rep(NA, times = length(un))
v9.pr2 <- data.frame(acc.tax = un, db = rep("pr2", times = length(un)), start = nav, stop = nav, asv.len = nav, asv = nav, stringsAsFactors = FALSE)
bloop <- c()
for (i in 1:length(un)) {
  fhits <- fsub[fsub$namer == un[i],]
  rhits <- rsub[rsub$namer == un[i],]
  strt <- fhits$fpos.end + 1
  fin <- rhits$rpos.strt - 1
  
  if (length(strt) == 1 && length(fin) == 1) {
    
    if (fin <= strt) {
      bloop <- append(bloop, i)
      v9.pr2[i, c("start","stop","asv.len")] <- c(strt, fin, length(strt:fin))
    } else {
      v9.pr2[i, c("start","stop","asv.len")] <- c(strt, fin, length(strt:fin))
      ss <- pr2v[un[i]]
      v9.pr2[i, "asv"] <- str_sub(ss, start = strt, end = fin)
    }
    
  } else if (length(strt) == 1 && length(fin) > 1) {
    # account for multiple reverse hits
    # use the amplicon who's size is closer to target
    eh <- which((fin - strt) == min(fin - strt))
    fin <- fin[eh]
    v9.pr2[i, c("start","stop","asv.len")] <- c(strt, fin, length(strt:fin))
    ss <- pr2v[un[i]]
    v9.pr2[i, "asv"] <- str_sub(ss, start = strt, end = fin)
    
  } else if (length(strt) > 1 && length(fin) == 1) {
    shit
    # account for multiple F hits
    # this never happens here...
  } else if (length(strt) > 1 && length(fin) > 1) {
    shyte
    # account for multiple F + R hits
    # this never happens here either...
  }
  
}

# remove bloop (where reverse hit further 5' than fwd) --> these are NA's so just rm NA's
v9.pr2 <- v9.pr2[which(!is.na(v9.pr2$asv)) , ]

### repeat above for silva-v9:
fsub <- silva.prime[[which(str_detect(names(silva.prime), paste0("fwd.", hv)))]]
fsub <- fsub[which(!is.na(fsub$namer)),]
rsub <- silva.prime[[which(str_detect(names(silva.prime), paste0("rev.", hv)))]]
rsub <- rsub[which(!is.na(rsub$namer)),]
silvav <- as.character(silva)

bb <- setdiff(fsub$namer, rsub$namer)
bloop <- fsub[fsub$namer %in% bb,]
strt <- bloop$fpos.end + 1
ss <- silvav[bloop$namer]
ll <- str_length(ss) - (strt+1)

v9.silva.noR <- data.frame(acc.tax = bloop$namer, db = rep("silva", times = length(bloop$namer)), start = strt, stop = str_length(ss), 
                         asv.len = str_length(ss) - (strt+1), asv = rep(NA, times = length(bloop$namer)), stringsAsFactors = FALSE)
for (i in 1:nrow(v9.silva.noR)) {
  v9.silva.noR[i, "asv"] <- str_sub(ss[i], start = v9.pr2.noR[i, "start"], end = length(ss))
}

bb <- intersect(fsub$namer, rsub$namer)
fsub <- fsub[fsub$namer %in% bb,]
rsub <- rsub[rsub$namer %in% bb,]
un <- unique(fsub$namer)
nav <- rep(NA, times = length(un))
v9.silva <- data.frame(acc.tax = un, db = rep("silva", times = length(un)), start = nav, stop = nav, asv.len = nav, asv = nav, stringsAsFactors = FALSE)
bloop <- c()
for (i in 1:length(un)) {
  fhits <- fsub[fsub$namer == un[i],]
  rhits <- rsub[rsub$namer == un[i],]
  strt <- fhits$fpos.end + 1
  fin <- rhits$rpos.strt - 1
  
  if (length(strt) == 1 && length(fin) == 1) {
    
    if (fin <= strt) {
      bloop <- append(bloop, i)
      v9.silva[i, c("start","stop","asv.len")] <- c(strt, fin, length(strt:fin))
    } else {
      v9.silva[i, c("start","stop","asv.len")] <- c(strt, fin, length(strt:fin))
      ss <- silvav[un[i]]
      v9.silva[i, "asv"] <- str_sub(ss, start = strt, end = fin)
    }
    
  } else if (length(strt) == 1 && length(fin) > 1) {
    fock
    # account for multiple reverse hits
    # this never happens here...
  } else if (length(strt) > 1 && length(fin) == 1) {
    fuck
    # account for multiple F hits
    # this never happens here...
  } else if (length(strt) > 1 && length(fin) > 1) {
    fack
    # account for multiple F + R hits
    # this never happens here...
  }
  
}
foock
# got all the v9 shit done and stored in 4 datasets
# cleaning should probably happen -- deal w/ 
# 1. multiple hitting primers in no R data set
# 2. weird ASV lengths 
# 3. ASV's containing degeneracies 
# 4. merging the 2 db's
# 5. check that hypervariable ASVs are still unique and trim to a common taxonomy if not (like LCA style)...

feck
# for now, just rm the ones w/ multiple primer hits on the same ref:
v9.pr2 <- v9.pr2[!(duplicated(v9.pr2$acc.tax) | duplicated(v9.pr2$acc.tax, fromLast = TRUE)), ]
v9.pr2.noR <- v9.pr2.noR[!(duplicated(v9.pr2.noR$acc.tax) | duplicated(v9.pr2.noR$acc.tax, fromLast = TRUE)), ]
v9.silva <- v9.silva[!(duplicated(v9.silva$acc.tax) | duplicated(v9.silva$acc.tax, fromLast = TRUE)), ]
v9.silva.noR <- v9.silva.noR[!(duplicated(v9.silva.noR$acc.tax) | duplicated(v9.silva.noR$acc.tax, fromLast = TRUE)), ]


# below shows that asv length has a bunch of NA's in it.. I think
# look at ASV lengths from the 4 data sets:
boxplot(c(v9.pr2$asv.len, v9.pr2.noR$asv.len, v9.silva$asv.len, v9.silva.noR$asv.len))
# remove ASVs longer than 180nt and shorter than 90nt (target ~120-130)
v9.pr2 <- v9.pr2[v9.pr2$asv.len > 90 & v9.pr2$asv.len < 180 ,]
v9.pr2.noR <- v9.pr2.noR[v9.pr2.noR$asv.len > 90 & v9.pr2.noR$asv.len < 180 ,]
v9.silva <- v9.silva[v9.silva$asv.len > 90 & v9.pr2$silva.len < 180 ,]
v9.silva.noR <- v9.silva.noR[v9.silva.noR$asv.len > 90 & v9.silva.noR$asv.len < 180 ,]
# check that it worked:
boxplot(c(v9.pr2$asv.len, v9.pr2.noR$asv.len, v9.silva$asv.len, v9.silva.noR$asv.len))

# here in writing
feck



# remove asv's that contain a degeneracy:
# v9.pr2[str_detect(v9.pr2$asv,"N") ,]

# this removes duplicated v9 ASVs, but I don't think I should do that entirely...?
# v9.pr2 <- v9.pr2[!(duplicated(v9.pr2$asv) | duplicated(v9.pr2$asv, fromLast = TRUE)), ]




# below doesn't really work... so try something else
# nav <- rep(NA, times = length(c(names(silva), names(pr2))))
# db <- c(rep("silva", times = length(names(silva))), rep("pr2",times = length(names(pr2))))
# v9all <- data.frame(acc.tax = c(names(silva), names(pr2)), db = db, n.fwdhit = nav, n.revhit = nav, asv = nav, asv.start = nav, asv.end = nav)
# v4all <- data.frame(acc.tax = c(names(silva), names(pr2)), db = db, n.fwdhit = nav, n.revhit = nav, asv = nav, asv.start = nav, asv.end = nav)
# v45all <- data.frame(acc.tax = c(names(silva), names(pr2)), db = db, n.fwdhit = nav, n.revhit = nav, asv = nav, asv.start = nav, asv.end = nav)

# for (i in 1:length(nav)) {
#   an <- v9all$acc.tax[i]
# 
#   eh <- lapply(pr2.prime, function(x) which(x$namer == an))
#   eh2 <- lapply(silva.prime, function(x) which(x$namer == an))
#   if (length(unlist(eh)) == 0 && length(unlist(eh2)) == 0) {
#     break
#   } else if (length(unlist(eh)) > 0) {
#     # manipulate pr2
#     eh <- unlist(lapply(eh, length))
#     eh <- names(eh)[eh > 0]
# 
#     if (str_detect(eh, "v9") && str_detect(eh, "v4") && str_detect(eh, "v4.5") ) {
#       # all 3 primers hit...
#       fsub <- pr2.prime[["fwd.pr2"]]
#       rsub <- pr2.prime[["rev.pr2"]]
#       feck
#     } else if (str_detect(eh, "v9")) {
#       # looks like this never happens
#       
#     } 
# 
#   } else if (length(unlist(eh2)) > 0) {
#     # manipulate silva
#   }
# }

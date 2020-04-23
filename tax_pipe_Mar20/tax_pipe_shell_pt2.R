# picking up where the first one left off (after mapping all tax tables onto the pr2 taxonomy)

# scrapping the whole thing and finding macro's much more inteligently, lol. Line 190ish
# ^nvm. wasn't accounting for NA's 

# Up next:
# 1. go back and bring in corrected lca/bayes-silva arrays
# 2. save prok results to .csvs
# 3. figure out what to do with macroeuk's...
# 3a. leaning -- if majority are non-protist, rm; if majority are protist or NA, keep til later...
# 4. also would be wise to save ASV seqs that are removed and/or questionable for each categorization
# 4a. this will make it easier to go back and analyze and/or remove them later

rm(list=setdiff(ls(), c("xx")))
setwd("~/Documents/R/amplicon_bioinformatics/tax_pipe_Mar20/")
source("~/Documents/R/amplicon_bioinformatics/package_deal/all_of_it.R")
library("stringr")
# xx <- readRDS(file = "all_mapped_taxtabs.rds")

orderByASV <- function(x) {
  ii <- base::sort(x$svN, index.return = TRUE)
  x <- x[ii$ix,]
}

xx <- lapply(xx, orderByASV)
nn <- names(xx)

# use silva tables to ID prok's (consistently and inconsistently across tables)
# silvas <- xx[str_which(nn,"silva")]
# silva.compz <- compare_by_tax_name(silvas[[1]], silvas[[2]], silvas[[3]], 
#                                    taxnames = c("Archaea", "Bacteria"),
#                                    tablenames = names(silvas), return.conflix = FALSE)
# arc.summary <- silva.compz[["Archaea"]][[1]]
# there are no conflicts for Archaea classifications:
# > arc.summary
# bayes.silva idtax.silva lca.silva Freq
# 1     Archaea     Archaea   Archaea  105
# 2        <NA>     Archaea   Archaea    5
# 3     Archaea        <NA>   Archaea    1
# 4        <NA>        <NA>   Archaea    4
# 6        <NA>     Archaea      <NA>    4
# 7     Archaea        <NA>      <NA>    7
# bac.summary <- silva.compz[["Bacteria"]][[1]]
# there are 20 total conflicts (ASVs assigned to something else - Eukaryota) for Bacteria
# > bac.summary
# bayes.silva idtax.silva lca.silva Freq
# 1     Bacteria    Bacteria  Bacteria 1729
# 3         <NA>    Bacteria  Bacteria    6
# 7     Bacteria        <NA>  Bacteria   44
# 8    Eukaryota        <NA>  Bacteria    1
# 9         <NA>        <NA>  Bacteria   32
# 13    Bacteria   Eukaryota Eukaryota    1
# 16    Bacteria        <NA> Eukaryota    2
# 19    Bacteria    Bacteria      <NA>   45
# 20   Eukaryota    Bacteria      <NA>    8
# 21        <NA>    Bacteria      <NA>    5
# 22    Bacteria   Eukaryota      <NA>    6
# 25    Bacteria        <NA>      <NA>  243

# bloop <- which((rowSums(bac.summary == "Bacteria", na.rm = TRUE) + rowSums(is.na(bac.summary))) < (ncol(bac.summary)-1))
# bac.i <- silva.compz[["Bacteria"]][[2]]
# bac.no.conflix.i <- unlist(bac.i[,-bloop]); bac.no.conflix.i <- bac.no.conflix.i[!is.na(bac.no.conflix.i)] # indices to remove as bacteria...
# arc.i <- silva.compz[["Archaea"]][[2]]
# arc.no.conflix.i <- unlist(arc.i[,-bloop]); arc.no.conflix.i <- arc.no.conflix.i[!is.na(arc.no.conflix.i)] # indices to remove as bacteria...
# 
# conflix.i <- unlist(bac.i[,bloop]); conflix.i <- conflix.i[!is.na(conflix.i)]
# # compute relative abundances of conflix:
# seqtab <- readRDS("seqtab_nochime_Mar20.rds")
# getme <- silvas[[1]]$ASV[conflix.i]
# 
# seqtab.r <- seqtab/rowSums(seqtab)
# conflix.ra <- seqtab.r[,getme]
# n0 <- rowSums(conflix.ra == 0) # number of 0's in each row
# conf.asvs.per.sample <- length(conflix.i) - n0
# conf.ra.per.sample <- rowSums(conflix.ra)
# # eh <- apply(conflix.ra, MARGIN = 2, FUN = max)
# # eh2 <- rowSums(conflix.ra) # sum of conflicting ASVs by sample...
# # # eh2[eh2 > 0]
# 
# prok.i.no.conflict <- c(bac.no.conflix.i,arc.no.conflix.i)
# getme <- silvas[[1]]$ASV[prok.i.no.conflict]
# # compute ASVs and reads removed x sample and x(sample x asv) by prok filter:
# proksub.noconflict <- seqtab.r[,getme]
# n0 <- rowSums(proksub.noconflict == 0) # number of 0's in each row
# rmproks.asvs.per.sample <- length(prok.i.no.conflict) - n0
# proks.ra.per.sample <- rowSums(proksub.noconflict)
# 
# # write the conflicting and conflicting + total prok results to .csv's for records
# ml <- max(length(silvas[[1]]$ASV[prok.i.no.conflict]),
#            length(silvas[[1]]$ASV[conflix.i]))
# rm.prok.ASVs <- data.frame(proks.no.conflict = c(silvas[[1]]$ASV[prok.i.no.conflict],rep(NA,times = ml-length(prok.i.no.conflict))),
#                            proks.w.conflict = c(silvas[[1]]$ASV[conflix.i],rep(NA,times = ml-length(conflix.i))))
# write.csv(rm.prok.ASVs, file = "tax_filtering_results/ASVs_rm_as_prok.csv")
# 
# if ( all(names(rmproks.asvs.per.sample) == names(proks.ra.per.sample)) & 
#      all(names(rmproks.asvs.per.sample) == names(conf.asvs.per.sample)) &
#      all(names(rmproks.asvs.per.sample) == names(conf.ra.per.sample)) ){
#   prok.summary <- cbind(rmproks.asvs.per.sample, conf.asvs.per.sample, 
#                         rmproks.asvs.per.sample + conf.asvs.per.sample,
#                         proks.ra.per.sample, conf.ra.per.sample,
#                         proks.ra.per.sample + conf.ra.per.sample)
#   colnames(prok.summary) <- c("rm.prok.asvs.no.conflict", "rm.prok.asvs.w.conflict", "total.rm.prok.asv",
#                               "rm.prok.rel.ab.no.conflict", "rm.prok.rel.ab.w.conflict", "total.rm.prok.rel.ab")
#   prok.summary <- as.data.frame(prok.summary, stringsAsFactors = FALSE)
# }
# write.csv(prok.summary, file = "tax_filtering_results/prok_per_sample_summary.csv")
# 
# # remove all potential proks
# rm.prok <- c(silvas[[1]]$ASV[prok.i.no.conflict],silvas[[1]]$ASV[conflix.i])
# # helper function to remove ASVs from each table in the list:
# rmmer <- function(x){
#   x <- x[!x$ASV %in% rm.prok,]
# }
# xx <- lapply(xx, FUN = rmmer)

# macro-euks:
macro.names <- c("Metazoa", "Fungi", "Streptophyta", "Rhodophyta", "Ulvophyceae", "Phaeophyceae")
macro.compz <- compare_by_tax_name(xx[[1]],xx[[2]],xx[[3]],xx[[4]],xx[[5]],xx[[6]],
                                   taxnames = macro.names,
                                   tablenames = names(xx),
                                   return.conflix = FALSE)
# write a loop to summarize
meta.summary <- macro.compz[["Metazoa"]][[1]]
fungi.summary <- macro.compz[["Fungi"]][[1]]
str.summary <- macro.compz[["Streptophyta"]][[1]]
rho.summary <- macro.compz[["Rhodophyta"]][[1]]
ulv.summary <- macro.compz[["Ulvophyceae"]][[1]]
pha.summary <- macro.compz[["Phaeophyceae"]][[1]]
rm.macro.maj <- c() # where macro were assigned as a majority (some protists still)
rm.macro.all <- c() # where macro was exclusively assigned (no possible protists, NA's allowed)
macro.minority <- c() # where majority of names are protist
get.propo <- function(x) {
  l <- length(which(x %in% macro.names))
  nna <- length(which(is.na(x)))
  nnamed <- (length(x)-1) - nna
  propo <- l / nnamed
  return(propo)
}
for (i in 1:length(macro.compz)) {
  ii <- macro.compz[[i]][[2]] # index array for this comp
  SUM <- macro.compz[[i]][[1]] # summary array for this comp
  yaboi <- apply(SUM, MARGIN = 1, FUN = get.propo)
  # names(yaboi) <- NULL
  # indices of yaboi > 0.5 & < 1 are cols of ii that need to be combined and added to rm.macro.maj
  # indices of yaboi == 1 are cols of ii that need to be combined and added to rm.macro.all
  # indices of yaboi < 0.5 are majority protist assigned -- keep them in but save for later
  i.maj <- which(yaboi > 0.5 & yaboi < 1)
  i.mall <- which(yaboi == 1)
  i.mino <- which(yaboi <= 0.5)
  # ASV indices in tax-tables for each category
  maj <- unlist(ii[,i.maj])
  maj <- maj[!is.na(maj)] # majority ASV (row) indices in taxtabs
  mall <- unlist(ii[,i.mall])
  mall <- mall[!is.na(mall)] # all ''
  mino <- unlist(ii[,i.mino])
  mino <- mino[!is.na(mino)] # min ''

  # # trying to troubleshoot but giving up...
  # if(any(maj %in% c(3626, 8991)) || any(mall %in% c(3626, 8991))){
  #   # if i == 2 these are in maj but not mall
  #   # not in rm.macro.all before proceeding
  #   # still not in rm.macro.all after proceeding
  #   # none of the other if's break either
  #
  #   # if i == 3 these are not in maj but in mall
  #   dam
  # }
  # if (any(c(2775,6581,6587,13980) %in% maj) || any(c(2775,6581,6587,13980) %in% mino)) {
  #   # if i == 1 these are in maj & not mino
  #   # also not in macro.minority if u continue...
  #   yikes
  # }
  # if (any(c(13006, 13202, 13966, 14178, 14349) %in% mino) || any(c(13006, 13202, 13966, 14178, 14349) %in% mall)) {
  #   ah
  # }

  rm.macro.maj <- unique(append(rm.macro.maj, maj))
  rm.macro.all <- unique(append(rm.macro.all, mall))
  macro.minority <- unique(append(macro.minority, mino))
}
# collect weird overlapping ones
weird <- c(intersect(rm.macro.maj,rm.macro.all),
           intersect(rm.macro.maj,macro.minority),
           intersect(rm.macro.all,macro.minority))
yy <- lapply(xx, function(x) x <- x[weird,])
y2 <- matrix(unlist(yy), nrow = length(weird), byrow = FALSE)
View(y2) # you can look at the weird ones here...
# after manual review, all the weird overlapping ones are majority [1:6] or all [7:end] non-protist 
macro.minority <- macro.minority[!macro.minority %in% weird]
rm.macro.maj <- rm.macro.maj[!rm.macro.maj %in% weird[7:length(weird)]]
rm.macro.all <- rm.macro.all[!rm.macro.all %in% weird[1:6]]

# this seems to work


# I thought this would be easier but it's not...

# macro.names <- c("Metazoa", "Fungi", "Streptophyta", "Rhodophyta", "Ulvophyceae", "Phaeophyceae")
# containMacro <- function(x) {
#   macrows <- apply(x, MARGIN = 1, function(z) any(z %in% macro.names))
#   ci <- apply(x[macrows,], MARGIN = 1, function(z) which(z))
#   return(list(macrows,ci))
# }
# # given a logical matrix, the following calcs proportion of TRUE elements in each row
# get.propo <- function(x) {
#   tt <- apply(x, MARGIN = 1, function(z) length(which(z)))
#   ff <- apply(x, MARGIN = 1, function(z) length(which(!z)))
#   propo <- tt/(tt+ff)
#   return(propo)
# }
# countNA <- function(x) {
#   
# }
# eh <- lapply(xx, containMacro)
# ci <- eh[[2]]
# eh
# ehm <- matrix(unlist(eh[[1]]), nrow = nrow(xx[[1]]), ncol = length(xx), byrow = FALSE)
# macro.propz <- get.propo(ehm)
# i.allmac <- which(macro.propz == 1)
# i.mostmac <- which(macro.propz > 0.5 & macro.propz < 1)
# i.somemac <- which(macro.propz > 0 & macro.propz <= 0.5)
# 
# # do this to check your work... (change the number in [3])
# # > lapply(xx, function(x) x[i.allmac[3],-c(1,2)])
# # > lapply(xx, function(x) x[i.mostmac[3],-c(1,2)])
# # > lapply(xx, function(x) x[i.somemac[3],-c(1,2)])
# 
# 
# prok.names <- c("Bacteria","Archaea")
# containProk <- function(x) {
#   macrows <- apply(x, MARGIN = 1, function(z) any(z %in% "prok.names"))
#   return(macrows)
# }
# eh <- lapply(xx, containProk)
# ehm <- matrix(unlist(eh), nrow = nrow(xx[[1]]), ncol = length(xx), byrow = FALSE)
# prok.propz <- get.propo(ehm)
# iprok <- which(prok.propz > 0)



# remove macro's ID'ed as all and/or as majority

feck
# note that it doesn't really make sense to do this until after removing non-protists...
# analyze the trait mapping results:
trts <- c("Chloroplast", "Plast_Origin", "Ingestion", "Cover", "Shape")
for (i in 1:length(xx)){
  tt <- xx[[i]]
  mr <- read.csv(paste0("unmappedtax_traitmap_results/",nn[i],"_mapout.csv"), stringsAsFactors = FALSE)
  for (j in 1:length(trts)){
    fout <- c(paste0("unmappedtax_traitmap_results/", nn[i], "_", trts[j], "_ASV_plot.pdf"), 
              paste0("unmappedtax_traitmap_results/", nn[i], "_", trts[j], "_abundance_plot1.pdf"),
              paste0("unmappedtax_traitmap_results/", nn[i], "_", trts[j], "_abundance_plot2.pdf"))
    analyze_traitmap_byTrait(map.result = mr, trait.name = trts[j], otu.table = seqtab, pltfilez = fout)
  }
}


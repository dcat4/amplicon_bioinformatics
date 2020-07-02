# uses the output from pt2 (the mapped, tax-filtered taxtabs) to:

# Kevin's algorithm isn't working properly... sv10002 is a good example. maybe an issue with tie-breaking
# issue appears to be at line 72ish - if the tax to favor is NA, it sets it to "na", which is the opposite of intended behavior

# to fix the consensus bugs, you should just run it for sv10002, return in the middle to see what the data looks like, and rectify from there
# or remove the tie-breaking and just weight tables differently

# traitmap analysis works but results are pretty terrible...

rm(list=ls())
setwd("~/Documents/R/amplicon_bioinformatics/tax_pipe_Mar20/")
source("~/Documents/R/amplicon_bioinformatics/package_deal/all_of_it.R")
xx <- readRDS(file = "all_mapped_taxtabs_protistOnly.rds")

seqtab <- readRDS("seqtab_nochime_Mar20.rds")
samsubs <- readRDS("sample_subsets.rds")

# check that ASVs have been removed from taxtab's, and remove the sames ones from seqtab
prok.asvs <- read.csv(file = "tax_filtering_results/ASVs_rm_as_prok.csv", stringsAsFactors = FALSE)
prok.asvs <- unlist(prok.asvs[, colnames(prok.asvs) %in% c("proks.no.conflict", "proks.w.conflict")], use.names = FALSE)
mac.asvs <- read.csv(file = "tax_filtering_results/ASVs_rm_as_macroEuk.csv", stringsAsFactors = FALSE)
mac.asvs <- unlist(mac.asvs[, colnames(mac.asvs) %in% c("rm.macro.all", "rm.macro.maj")], use.names = FALSE)
rm.asvs <- c(prok.asvs, mac.asvs)
rm.asvs <- rm.asvs[!is.na(rm.asvs)]

isitthere <- function(x) {
  eh <- any(rm.asvs %in% x)
}
hmm <- any(lapply(xx, isitthere))
if (hmm) {
  stop("asvs remain in taxtables. revisit pt 2")
}

seqtab <- seqtab[, !(colnames(seqtab) %in% rm.asvs)]

# can add optional sample subsets here...


## tax table comparisons:
# pairwise across algorithms:
rn <- c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species")
tblnam <- c("bayes-pr2", "idtax-pr2", "lca-pr2")
pr2.3way.byrank <- compare_byRank_3way(xx[[tblnam[1]]],xx[[tblnam[2]]],xx[[tblnam[3]]],
                                       pltfilez = "none",
                                       tablenames = tblnam,
                                       ranknamez = rn)
tblnam <- c("bayes-silva", "idtax-silva", "lca-silva")
silva.3way.byrank <- compare_byRank_3way(xx[[tblnam[1]]],xx[[tblnam[2]]],xx[[tblnam[3]]],
                                       pltfilez = "none",
                                       tablenames = tblnam,
                                       ranknamez = rn)
# cowplot and save 3-way plots:
library("cowplot")
yl <- c(0,1)
lowW <- 0.5 # reducing widths by this proportion for broader summary panels
p.3way.broad <- plot_grid(
  pr2.3way.byrank[[4]] + coord_cartesian(ylim = yl), #+ theme(legend.position = "none") #+ ggtitle("pr2 w/ 3 algorithms") ,
  # pr2.3way.byrank[[5]] + ggtitle("") + coord_cartesian(ylim = yl),
  silva.3way.byrank[[4]] + coord_cartesian(ylim = yl), #+ theme(legend.position = "none") + ggtitle("silva w/ 3 algorithms") ,
  # silva.3way.byrank[[5]] + ggtitle("") + coord_cartesian(ylim = yl),
  align = 'hv',
  labels = c("A.", "B."),
  # rel_widths = c(lowW, 1, lowW, 1),
  axis = 'l',
  hjust=-1,
  nrow = 2
)
ggsave("mapped_taxtab_comp_results/broad_bothDBs_3way_byDB_byRank.pdf", p.3way.broad, width = 8, height = 11.5, units = "in", device="pdf")

p.3way.deet <- plot_grid(
  pr2.3way.byrank[[5]] + coord_cartesian(ylim = yl),
  silva.3way.byrank[[5]] + coord_cartesian(ylim = yl),
  align = 'hv',
  labels = c("A.", "B."),
  axis = 'l',
  hjust=-1,
  nrow = 2
)
ggsave("mapped_taxtab_comp_results/deet_bothDBs_3way_byDB_byRank.pdf", p.3way.deet, width = 11, height = 13, units = "in", device="pdf")

# do by 2-ways by algorithm across each database
tblnam <- c("bayes-pr2", "bayes-silva")
bayes.2way.byRank <- compare_byRank_2way(xx[[tblnam[1]]],xx[[tblnam[2]]],
                                         pltfilez = "none",
                                         tablenames = tblnam,
                                         ranknamez = rn)
# bayes.2way.byAss <- compare_assignments_2way(xx[[tblnam[1]]],xx[[tblnam[2]]],
#                                          pltfilez = "none",
#                                          tablenames = tblnam,
#                                          ranknamez = rn)
tblnam <- c("idtax-pr2", "idtax-silva")
idtax.2way.byRank <- compare_byRank_2way(xx[[tblnam[1]]],xx[[tblnam[2]]],
                                         pltfilez = "none",
                                         tablenames = tblnam,
                                         ranknamez = rn)
# idtax.2way.byAss <- compare_assignments_2way(xx[[tblnam[1]]],xx[[tblnam[2]]],
#                                              pltfilez = "none",
#                                              tablenames = tblnam,
#                                              ranknamez = rn)
tblnam <- c("lca-pr2", "lca-silva")
lca.2way.byRank <- compare_byRank_2way(xx[[tblnam[1]]],xx[[tblnam[2]]],
                                         pltfilez = "none",
                                         tablenames = tblnam,
                                         ranknamez = rn)
# lca.2way.byAss <- compare_assignments_2way(xx[[tblnam[1]]],xx[[tblnam[2]]],
#                                              pltfilez = "none",
#                                              tablenames = tblnam,
#                                              ranknamez = rn)
# cowplot dem 2-way fools:
# p.2way <- plot_grid(
#   bayes.2way.byRank[[4]] + ggtitle("by Rank\nbayes") + coord_cartesian(ylim = yl),
#   # bayes.2way.byAss[[3]] + ggtitle("by Assignment\nbayes") + coord_cartesian(ylim = yl),
#   idtax.2way.byRank[[4]] + ggtitle("idtax") + coord_cartesian(ylim = yl),
#   # idtax.2way.byAss[[3]] + ggtitle("idtax") + coord_cartesian(ylim = yl),
#   lca.2way.byRank[[4]] + ggtitle("LCA") + coord_cartesian(ylim = yl),
#   # lca.2way.byAss[[3]] + ggtitle("LCA") + coord_cartesian(ylim = yl),
#   align = 'hv',
#   # labels = c("A.", "B.", "C.","D.","E.","F."),
#   axis = 'l',
#   hjust=-1,
#   nrow = 3
# )
# ggsave("mapped_taxtab_comp_results/byAlgorithm_2way.pdf", p.2way, width = 13, height = 18, units = "in", device="pdf")

# a check for the frankenstein algorithm:
# yy <- yy <- c(xx, xx[1])
# names(yy) <- c(names(xx), "ensemble")
# check4frankenstein(yy)
# the above should print "no frankenstein assignments detected" and it does so u good bro

# the below seems to work great/// frankenstein checking is slow but fine.

# create 2 ensembles (NA or not) and compare with each individual table:
all.c <- consensus_tax_mostCom2(xx[[1]],xx[[2]],xx[[3]],xx[[4]],xx[[5]],xx[[6]],
                               tablenames=names(xx), 
                               ranknamez=colnames(xx[[1]])[3:length(colnames(xx[[1]]))],
                               tiebreakz=list(c("idtax-pr2", NA), c("idtax-silva", NA),c("bayes-pr2", NA), c("bayes-silva", NA),c("lca-pr2", NA)), 
                               count.na=TRUE, 
                               trueMajority=FALSE,
                               weights=c(1,1,1,1,1,1))
yy <- c(xx, list(all.c))
names(yy) <- c(names(xx), "ensemble")
check4frankenstein(yy)
feck

## add sample subsetting routine here...

## traitmap analysis for unmapped tax-tables (removes macro's and proks)
# works but useless without sample subset context...
# i only had ASVs w/ the LCA arrays during mapping, but all taxtabs had ASVs in same order
# so pop the ASVs out from LCA and put them on all the others
tt <- read.csv(paste0("unmappedtax_traitmap_results/lca-pr2_mapout.csv"), stringsAsFactors = FALSE)
map.sv <- tt[,c("svN","ASV")]
rm.me <- which(map.sv$ASV %in% rm.asvs)
## read in pre-mapped traitmapping results and work them up
nn <- names(xx)
# note that it doesn't really make sense to do this until after removing non-protists...
# analyze the trait mapping results:
# trts <- c("Chloroplast", "Plast_Origin", "Ingestion", "Cover", "Shape")
trts <- c("Chloroplast")
for (k in 1:ncol(samsubs)){
  ss <- samsubs[,1]
  ss <- ss[!is.na(ss)]
  sgrp <- colnames(samsubs)[k]
  for (i in 1:length(nn)){
    ot <- seqtab[rownames(seqtab) %in% ss , ]
    mr <- read.csv(paste0("unmappedtax_traitmap_results/",nn[i],"_mapout.csv"), stringsAsFactors = FALSE)
    mr <- cbind(map.sv, mr)
    mr <- mr[-rm.me,]
    nothere <- which(colSums(ot) == 0) # ASVs not detected in this sample group
    nothere <- colnames(ot)[nothere]
    mr <- mr[!mr$ASV %in% nothere ,]
    ot <- ot[, !colnames(ot) %in% nothere]
    for (j in 1:length(trts)){
      fout <- c(paste0("unmappedtax_traitmap_results/", sgrp, "/", nn[i], "_", trts[j], "_ASV_plot.pdf"), 
                paste0("unmappedtax_traitmap_results/", sgrp, "/", nn[i], "_", trts[j], "_abundance_plot1.pdf"),
                paste0("unmappedtax_traitmap_results/", sgrp, "/", nn[i], "_", trts[j], "_abundance_plot2.pdf"))
      analyze_traitmap_byTrait(map.result = mr, trait.name = trts[j], otu.table = ot, pltfilez = fout)
    }
  }
}

# mapping tax-mapped and filtered tables to trait db
# mm <- readRDS(file = "~/Documents/R/amplicon_bioinformatics/taxonomy_pipeline/tax_table_mapping/Ramond_traitdb_clean.RDS")
# dm <- c("Eukaryota", "Archaea", "Bacteria", "Metazoa", "Fungi", "Streptophyta", "Rhodophyta", "Ulvophyceae", "Phaeophyceae",
#         "Alveolata","Opisthokonta","Archaeplastida","Excavata","Rhizaria","Stramenopiles",
#         "Hacrobia","Amoebozoa","Apusozoa","Eukaryota_X","Protalveolata","Terrabacteria")
# for (i in 1:length(mappedtt)) {
#   tt <- mappedtt[[i]]
#   fout <- c(paste0("mappedtax_traitmap_results/",nn[i],"_mapout.csv"), 
#             paste0("mappedtax_traitmap_results/",nn[i],"_namesNotMapped.csv"))
#   traitmapper_Ramond(tt, mm, dont.map = dm, filezout = fout)
# }
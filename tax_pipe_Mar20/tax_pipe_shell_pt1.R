# UNDER CONSTRUCTION

# up next:
# 1. add pre-processing to rm prok's, small ASVs, metazoans --> this doesn't really work until you decide bootstrap values...
# 1a. propagate find_asvs_by_name to taxonomy_pipeline
# 2. pick and apply boot-strap thresholds -- trying defaults for now...
# 3. prelim trait mapping and analysis

# this is a shell script that executes the functions I (+Kevin +Connie) written for my ensemble taxonomy pipeline
# also using it as an outline to track where I'm at from start to finish..
# more of just an outline right now...

#### Starting point:
# Starts with initial taxonomic assignments output by DADA2/idtaxa/BLAST-MEGAN-LCA, boot-strapping thresholds not set yet, with arrays/files as follows:
# 1. bayesian classifier vs. pr2 --> /initial_tax_tabs/LCA_pr2_rawdf_Mar20.rds
# 2. bayesian classifier vs. silva --> /initial_tax_tabs/LCA_silva_rawdf_Mar20.rds
# 3. idtaxa vs. silva --> /initial_tax_tabs/idtax_silva_0boot_Mar20.rds
# 4. idtaxa vs. pr2 --> /initial_tax_tabs/idtax_pr2_0boot_Mar20.rds
# 5. MEGAN_LCA vs. Silva --> /initial_tax_tabs/bayes_silva_0boot_Mar20.rds
# 6. MEGAN_LCA vs. pr2 --> /initial_tax_tabs/bayes_pr2_0boot_Mar20.rds

# a collection of R scripts and command line calls was used to generate the above, see: 
# 1. bayesian classifier vs. pr2 --> assign_bayesTax_0boot.R, assignment_shell.R	
# 2. bayesian classifier vs. silva --> assign_bayesTax_0boot.R, assignment_shell.R	
# 3. idtaxa vs. silva --> idtaxa_vs_silva_0boot.R, assignment_shell.R	
# 4. idtaxa vs. pr2 --> idtaxa_vs_pr2_0boot.R, assignment_shell.R	
# 5. MEGAN_LCA vs. Silva --> get_seqs4blast.R, bring_lca2R.R, files in /blaster
# 6. MEGAN_LCA vs. pr2 [in 2 chunks] --> gethalf_seqs4blast.R, bring_lca2R.R, files in /blaster

#### Step 0: formatting taxonomy tables for analysis

# clear your workspace, setwd, read in your fcns:
rm(list=ls())
setwd("~/Documents/R/amplicon_bioinformatics/tax_pipe_Mar20/")
source("~/Documents/R/amplicon_bioinformatics/package_deal/all_of_it.R")

# read in your initial tax tabs:
bayes.pr2 <- readRDS("initial_tax_tabs/bayes_pr2_0boot_Mar20.rds")
idtax.pr2 <- readRDS("initial_tax_tabs/idtax_pr2_0boot_Mar20.rds")
bayes.silva <- readRDS("initial_tax_tabs/bayes_silva_0boot_Mar20.rds")
idtax.silva <- readRDS("initial_tax_tabs/idtax_silva_0boot_Mar20.rds")
lca.pr2 <- readRDS("initial_tax_tabs/LCA_pr2_rawdf_Mar20.rds")
lca.silva <- readRDS("initial_tax_tabs/LCA_silva_rawdf_Mar20.rds")

# here's the rubric for aligning ASV numbers and sequences across datasets:
library("DECIPHER")
rubber <- readDNAStringSet("blaster/allASVs4blast.fasta")

# convert idtax and bayes to dataframes - keep boot = 0  and return conf for now...
# note that for idtaxa, order of ASVs must be retrieved from dada2 seqtab - rubric represents that order
# and adds ASV names and sequences to output dataframe:
xx <- bayestax2df(bayes.pr2, boot = 0, rubric = rubber, return.conf = TRUE)
bayes.pr2 <- xx[[1]]
bayes.pr2.conf <- xx[[2]]

xx <- bayestax2df(bayes.silva, boot = 0, rubric = rubber, return.conf = TRUE)
bayes.silva <- xx[[1]]
bayes.silva.conf <- xx[[2]]

xx <- idtax2df_pr2(idtax.pr2, boot = 0, rubric = rubber, return.conf = TRUE)
idtax.pr2 <- xx[[1]]
idtax.pr2.conf <- xx[[2]]

xx <- idtax2df_silva(idtax.silva, boot = 0, rubric = rubber, return.conf = TRUE)
idtax.silva <- xx[[1]]
idtax.silva.conf <- xx[[2]]

# lca.pr2 + lca.silva were already formatted when converting from .csv, so now we're good with formatting

# clear out some unnecessary stuff:
rm("rubber", "xx")

# re-order your 6 tables by sorting according to ASV:
ii <- base::sort(bayes.pr2$ASV, index.return = TRUE)
bayes.pr2 <- bayes.pr2[ii$ix,]
bayes.pr2.conf <- bayes.pr2.conf[ii$ix,]

ii <- base::sort(idtax.pr2$ASV, index.return = TRUE)
idtax.pr2 <- idtax.pr2[ii$ix,]
idtax.pr2.conf <- idtax.pr2.conf[ii$ix,]

ii <- base::sort(lca.pr2$ASV, index.return = TRUE)
lca.pr2 <- lca.pr2[ii$ix,]
# no confidence for lca tax assignments...

ii <- base::sort(bayes.silva$ASV, index.return = TRUE)
bayes.silva <- bayes.silva[ii$ix,]
bayes.silva.conf <- bayes.silva.conf[ii$ix,]

ii <- base::sort(idtax.silva$ASV, index.return = TRUE)
idtax.silva <- idtax.silva[ii$ix,]
idtax.silva.conf <- idtax.silva.conf[ii$ix,]

ii <- base::sort(lca.silva$ASV, index.return = TRUE)
lca.silva <- lca.silva[ii$ix,]
# no confidence for lca tax assignments...

# check that they're all in the same order
identical(bayes.pr2$ASV, bayes.silva$ASV)
identical(bayes.pr2$ASV, idtax.pr2$ASV)
identical(idtax.pr2$ASV, idtax.silva$ASV)  
identical(idtax.silva$ASV, lca.pr2$ASV)
identical(lca.pr2$ASV, lca.silva$ASV)

# another peak b/c i'm ocd:
head(bayes.pr2)
head(bayes.silva)
head(idtax.pr2)
head(idtax.silva)
head(lca.pr2)
head(lca.silva)
# Noice.

#### Step 1: boot threshold sensitivity analysis and preliminary mapping/comparisons

### 1a. bootstrap threshold comparison
## Here loop thru a gradient of bootstrapping thresholds and compare:
# resolution of assignments across the data set
# pairwise comparisons of =thresholds for idtaxa vs. bayes of each ref db
bayes.silva <- bayes.silva[, -which(colnames(bayes.silva) %in% c("svN", "ASV"))]
bayes.pr2 <- bayes.pr2[, -which(colnames(bayes.pr2) %in% c("svN", "ASV"))]
idtax.silva <- idtax.silva[, -which(colnames(idtax.silva) %in% c("svN", "ASV"))]
idtax.pr2 <- idtax.pr2[, -which(colnames(idtax.pr2) %in% c("svN", "ASV"))]
# vector of boot thresholds u want to fux wit
bootvec <- seq(from = 40, to = 80, by = 10)
bootvec.str <- vector(mode = "character")
# initialize lists for various boot thresholds applied to the 4 tables...
bayes.pr2.list <- rep(list(bayes.pr2), length(bootvec))
idtax.pr2.list <- rep(list(idtax.pr2), length(bootvec))
bayes.silva.list <- rep(list(bayes.silva), length(bootvec))
idtax.silva.list <- rep(list(idtax.silva), length(bootvec))
# 
# # for storing ggplots
# plot.list.pr2 <- list()
# plot.list.silva <- list()
# # this loop NAs-out assignments based on thresholds in boot. plots 2-way comparisons of pr2 and silva tax tabs at each boot thresholds
# for (i in 1:length(bootvec)) {
#   
#   bootvec.str <- append(bootvec.str, paste0("boot = ",toString(bootvec[i]),"%")) # for plot titles below...
#   
#   x1 <- bayes.pr2.list[[i]]
#   x1[bayes.pr2.conf < bootvec[i]] <- NA
#   bayes.pr2.list[[i]] <- x1
#   
#   x2 <- idtax.pr2.list[[i]]
#   x2[idtax.pr2.conf < bootvec[i]] <- NA
#   idtax.pr2.list[[i]] <- x2
#   
#   x3 <- bayes.silva.list[[i]]
#   x3[bayes.silva.conf < bootvec[i]] <- NA
#   bayes.silva.list[[i]] <- x3
#   
#   x4 <- idtax.silva.list[[i]]
#   x4[idtax.silva.conf < bootvec[i]] <- NA
#   idtax.silva.list[[i]] <- x4
# }
# 
# # will stitch your plots together (with cowplot) and save
# library("cowplot")
# yl <- c(0,1) # y-limits
# 
# # first look at defaults/recommended boots for each algorithm (idtax = 60, bayes = 80)
# # pr2:
# b80.pr2 <- bayes.pr2.list[[which(bootvec == 80)]]
# i60.pr2 <- idtax.pr2.list[[which(bootvec == 60)]]
# p1 <- compare_assignments_2way(b80.pr2, i60.pr2,
#                                pltfilez = "none",
#                                tablenames = c("bayes-pr2","idtax-pr2"),
#                                ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"))
# p2 <- compare_byRank_2way(b80.pr2, i60.pr2,
#                           pltfilez = "none",
#                           tablenames = c("bayes-pr2","idtax-pr2"),
#                           ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"))
# # silva:
# b80.silva <- bayes.silva.list[[which(bootvec == 80)]]
# i60.silva <- idtax.silva.list[[which(bootvec == 60)]]
# p3 <- compare_assignments_2way(b80.silva, i60.silva,
#                               pltfilez = "none",
#                               tablenames = c("bayes-silva","idtax-silva"),
#                               ranknamez = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))
# p4 <- compare_byRank_2way(b80.silva, i60.silva,
#                           pltfilez = "none",
#                           tablenames = c("bayes-silva","idtax-silva"),
#                           ranknamez = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))
# p.2way.default.bothdb <- plot_grid(
#   p1[[3]] + ggtitle("bayes-pr2, boot=80% vs. idtax-pr2, boot=60%") + coord_cartesian(ylim = yl),
#   p2[[4]] + coord_cartesian(ylim = yl),
#   p3[[3]] + ggtitle("bayes-silva, boot=80% vs. idtax-silva, boot=60%") + coord_cartesian(ylim = yl),
#   p4[[4]] + coord_cartesian(ylim = yl),
#   align = 'hv',
#   labels = c("A.", "B.", "C.","D."),
#   axis = 'l',
#   hjust=-1,
#   nrow = 2
# )
# ggsave("bayesVidtax_boot_threshold/bothdb_default_boots.pdf", p.2way.default.bothdb, width = 13, height = 12, units = "in", device="pdf")
# 
# # according to idtaxa creators, 50% is also "high confidence" so try that too...
# i50.pr2 <- idtax.pr2.list[[which(bootvec == 50)]]
# p5 <- compare_assignments_2way(b80.pr2, i50.pr2,
#                                pltfilez = "none",
#                                tablenames = c("bayes-pr2","idtax-pr2"),
#                                ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"))
# p6 <- compare_byRank_2way(b80.pr2, i50.pr2,
#                           pltfilez = "none",
#                           tablenames = c("bayes-pr2","idtax-pr2"),
#                           ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"))
# # silva:
# i50.silva <- idtax.silva.list[[which(bootvec == 50)]]
# p7 <- compare_assignments_2way(b80.silva, i50.silva,
#                                pltfilez = "none",
#                                tablenames = c("bayes-silva","idtax-silva"),
#                                ranknamez = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))
# p8 <- compare_byRank_2way(b80.silva, i50.silva,
#                           pltfilez = "none",
#                           tablenames = c("bayes-silva","idtax-silva"),
#                           ranknamez = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))
# p.2way.default.bothdb2 <- plot_grid(
#   p5[[3]] + ggtitle("bayes-pr2, boot=80% vs. idtax-pr2, boot=50%") + coord_cartesian(ylim = yl),
#   p6[[4]] + coord_cartesian(ylim = yl),
#   p7[[3]] + ggtitle("bayes-pr2, boot=80% vs. idtax-pr2, boot=50%") + coord_cartesian(ylim = yl),
#   p8[[4]] + coord_cartesian(ylim = yl),
#   align = 'hv',
#   labels = c("A.", "B.", "C.","D."),
#   axis = 'l',
#   hjust=-1,
#   nrow = 2
# )
# ggsave("bayesVidtax_boot_threshold/bothdb_default_boots2.pdf", p.2way.default.bothdb2, width = 13, height = 12, units = "in", device="pdf")
# 
# ### now adjust bayes-pr2 down iteratively, holding idtax at 50%
# # pr2:
# b70.pr2 <- bayes.pr2.list[[which(bootvec == 70)]]
# p1 <- compare_assignments_2way(b70.pr2, i50.pr2,
#                                pltfilez = "none",
#                                tablenames = c("bayes-pr2","idtax-pr2"),
#                                ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"))
# p2 <- compare_byRank_2way(b70.pr2, i50.pr2,
#                           pltfilez = "none",
#                           tablenames = c("bayes-pr2","idtax-pr2"),
#                           ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"))
# # silva:
# b70.silva <- bayes.silva.list[[which(bootvec == 70)]]
# p3 <- compare_assignments_2way(b70.silva, i50.silva,
#                                pltfilez = "none",
#                                tablenames = c("bayes-silva","idtax-silva"),
#                                ranknamez = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))
# p4 <- compare_byRank_2way(b70.silva, i50.silva,
#                           pltfilez = "none",
#                           tablenames = c("bayes-silva","idtax-silva"),
#                           ranknamez = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))
# p.2way.default.bothdb <- plot_grid(
#   p1[[3]] + ggtitle("bayes-pr2, boot=70% vs. idtax-pr2, boot=50%") + coord_cartesian(ylim = yl),
#   p2[[4]] + coord_cartesian(ylim = yl),
#   p3[[3]] + ggtitle("bayes-silva, boot=70% vs. idtax-silva, boot=50%") + coord_cartesian(ylim = yl),
#   p4[[4]] + coord_cartesian(ylim = yl),
#   align = 'hv',
#   labels = c("A.", "B.", "C.","D."),
#   axis = 'l',
#   hjust=-1,
#   nrow = 2
# )
# ggsave("bayesVidtax_boot_threshold/bothdb_default-10_boots.pdf", p.2way.default.bothdb, width = 13, height = 12, units = "in", device="pdf")
# 
# # according to idtaxa creators, 50% is also "high confidence" so try that too...
# b60.pr2 <- bayes.pr2.list[[which(bootvec == 60)]]
# p5 <- compare_assignments_2way(b60.pr2, i50.pr2,
#                                pltfilez = "none",
#                                tablenames = c("bayes-pr2","idtax-pr2"),
#                                ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"))
# p6 <- compare_byRank_2way(b60.pr2, i50.pr2,
#                           pltfilez = "none",
#                           tablenames = c("bayes-pr2","idtax-pr2"),
#                           ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"))
# # silva:
# b60.silva <- bayes.silva.list[[which(bootvec == 60)]]
# p7 <- compare_assignments_2way(b60.silva, i50.silva,
#                                pltfilez = "none",
#                                tablenames = c("bayes-silva","idtax-silva"),
#                                ranknamez = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))
# p8 <- compare_byRank_2way(b60.silva, i50.silva,
#                           pltfilez = "none",
#                           tablenames = c("bayes-silva","idtax-silva"),
#                           ranknamez = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))
# p.2way.default.bothdb2 <- plot_grid(
#   p5[[3]] + ggtitle("bayes-pr2, boot=60% vs. idtax-pr2, boot=50%") + coord_cartesian(ylim = yl),
#   p6[[4]] + coord_cartesian(ylim = yl),
#   p7[[3]] + ggtitle("bayes-pr2, boot=60% vs. idtax-pr2, boot=50%") + coord_cartesian(ylim = yl),
#   p8[[4]] + coord_cartesian(ylim = yl),
#   align = 'hv',
#   labels = c("A.", "B.", "C.","D."),
#   axis = 'l',
#   hjust=-1,
#   nrow = 2
# )
# ggsave("bayesVidtax_boot_threshold/bothdb_default-10or20_boots.pdf", p.2way.default.bothdb2, width = 13, height = 12, units = "in", device="pdf")
# 

# clear out the shit
rm(list=setdiff(ls(), c("bayes.pr2","bayes.pr2.conf","bayes.silva","bayes.silva.conf",
                        "idtax.pr2","idtax.pr2.conf","idtax.silva","idtax.silva.conf",
                        "lca.pr2", "lca.silva", "bootvec", "bootvec.str")))
# re-source your fcns:
source("~/Documents/R/amplicon_bioinformatics/package_deal/all_of_it.R")

# boot threshold implementation
bboot <- 60
iboot <- 50
bayes.pr2[bayes.pr2.conf < bboot] <- NA
bayes.silva[bayes.silva.conf < bboot] <- NA
idtax.pr2[idtax.pr2.conf < iboot] <- NA
idtax.silva[idtax.silva.conf < iboot] <- NA

## trying to write a new function to "compare_by_tax_name"
dummy <- compare_by_tax_name(bayes.silva, idtax.silva, lca.silva,
                             taxnames = c("Bacteria","Archaea"),
                             tablenames = c("bayes-silva","idtax-silva","lca-silva"),
                             return.conflix = FALSE)
dummy2 <- compare_by_tax_name(bayes.silva, idtax.silva, lca.silva,
                              taxnames = c("Bacteria","Archaea"),
                              tablenames = c("bayes-silva","idtax-silva","lca-silva"),
                              return.conflix = TRUE)
## this Bacteria thang all checks out...can just use dummy2 [maybe do one more check for metazoans]
# bac <- dummy$Bacteria
# bac.summ <- bac[[1]]; bac.i <- bac[[2]]
# # rows (cols) of bac.sum (bac.i) where there are conflicts (disagreeing names rather than target name+NA's)
# conflix.i <- which((rowSums(bac.summ == "Bacteria", na.rm = TRUE) + rowSums(is.na(bac.summ))) < (ncol(bac.summ)-1))
# rm.bac <- unlist(bac.i[,-conflix.i]); rm.bac <- rm.bac[!is.na(rm.bac)] # indices to remove as bacteria...
# checkme <- unlist(bac.i[,conflix.i]); checkme <- checkme[!is.na(checkme)]
# length(checkme) == sum(bac.summ$Freq[conflix.i])
# dummy2$Bacteria[[1]]
# dummy2$Bacteria[[2]]
# 
# arc <- dummy$Archaea
# arc.summ <- arc[[1]]; arc.i <- arc[[2]]
# # rows (cols) of bac.sum (bac.i) where there are conflicts (disagreeing names rather than target name+NA's)
# conflix.i <- which((rowSums(arc.summ == "Archaea", na.rm = TRUE) + rowSums(is.na(arc.summ))) < (ncol(arc.summ)-1))
# rm.arc <- unlist(arc.i[,-conflix.i]); rm.arc <- rm.arc[!is.na(rm.arc)] # indices to remove as bacteria...
# checkme <- unlist(arc.i[,conflix.i]); checkme <- checkme[!is.na(checkme)]
# length(checkme) == sum(arc.summ$Freq[conflix.i])

# attempting to remove proks:

#### Step 2: Mapping all tax tables to a common taxonomic nomenclature
# Kevin's stuff goes here

#### Step 3: Re-mapping individual tax-mapped tax tables onto trait database
# make sure boot-strapping thresholds are honored here...
# re-format as dataframes...
# remove non-protists...
# then write a shell of my trait mapping and analysis functions

#### Step 4: Comparisons of mapped taxonomy tables
# should be straightforward...

#### Step 5: ensemble taxonomy generation
# should be straightforward

#### Step 6: trait mapping of ensemble taxonomies
# should be straightforward


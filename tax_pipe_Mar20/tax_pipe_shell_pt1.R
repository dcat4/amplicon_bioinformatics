# UNDER CONSTRUCTION

# up next:
# 1. add pre-processing to rm prok's, small ASVs, metazoans --> this doesn't really work until you decide bootstrap values...
# 1a. propagate find_asvs_by_name to taxonomy_pipeline
# 2. pick and apply boot-strap thresholds
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

#### ADD PREPROCESSING HERE!!!
# library("stringr")
# # remove ASVs based on seq length (target = 120-130, trim < 90 & > 180):
# sl <- str_length(bayes.pr2$ASV)
# rmme <- which(sl < 90 | sl > 180)
# rmASV.seqlength <- bayes.pr2$ASV[rmme] # store it to summarize pre-processing later
# bayes.pr2 <- bayes.pr2[-rmme,]
# bayes.pr2.conf <- bayes.pr2.conf[-rmme,]
# bayes.silva <- bayes.silva[-rmme,]
# bayes.silva.conf <- bayes.silva.conf[-rmme,]
# idtax.pr2 <- idtax.pr2[-rmme,]
# idtax.pr2.conf <- idtax.pr2.conf[-rmme,]
# idtax.silva <- idtax.silva[-rmme,]
# idtax.silva.conf <- idtax.silva.conf[-rmme,]
# lca.pr2 <- lca.pr2[-rmme,]
# lca.silva <- lca.silva[-rmme,]
# 
# n2f.prok <- c("Bacteria", "Archaea")
# i.prok.all <- find_asvs_by_name(bayes.silva, idtax.silva, lca.silva, 
#                             names2find = n2f.prok, return.byTable = FALSE)
# rmme <- unlist(i.prok.all, use.names=FALSE)
# rmme <- rmme[!is.na(rmme)]
# rnASV.prok <- bayes.pr2$ASV[rmme]
# bayes.pr2 <- bayes.pr2[-rmme,]
# bayes.pr2.conf <- bayes.pr2.conf[-rmme,]
# bayes.silva <- bayes.silva[-rmme,]
# bayes.silva.conf <- bayes.silva.conf[-rmme,]
# idtax.pr2 <- idtax.pr2[-rmme,]
# idtax.pr2.conf <- idtax.pr2.conf[-rmme,]
# idtax.silva <- idtax.silva[-rmme,]
# idtax.silva.conf <- idtax.silva.conf[-rmme,]
# lca.pr2 <- lca.pr2[-rmme,]
# lca.silva <- lca.silva[-rmme,]
# 
# i.prok.lca <- find_asvs_by_name(bayes.silva, names2find = n2f.prok, return.byTable = FALSE)
# i.prok.bayes <- find_asvs_by_name(bayes.silva, names2find = n2f.prok, return.byTable = FALSE)
# i.prok.idtax <- find_asvs_by_name(idtax.silva, names2find = n2f.prok, return.byTable = FALSE)
# 
# feck

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

# for storing ggplots
plot.list.pr2 <- list()
plot.list.silva <- list()
# this loop NAs-out assignments based on thresholds in boot. plots 2-way comparisons of pr2 and silva tax tabs at each boot thresholds
for (i in 1:length(bootvec)) {
  
  bootvec.str <- append(bootvec.str, paste0("boot = ",toString(bootvec[i]),"%")) # for plot titles below...
  
  x1 <- bayes.pr2.list[[i]]
  x1[bayes.pr2.conf < bootvec[i]] <- NA
  bayes.pr2.list[[i]] <- x1
  
  x2 <- idtax.pr2.list[[i]]
  x2[idtax.pr2.conf < bootvec[i]] <- NA
  idtax.pr2.list[[i]] <- x2
  # 2-way comparison of pr2 tables
  r2way.pr2 <- compare_byRank_2way(x1, x2,
                                  pltfilez = "none",
                                  tablenames = c("bayes-pr2", "idtax-pr2"), 
                                  ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"))
  a2way.pr2 <- compare_assignments_2way(x1, x2, 
                                        pltfilez = "none",
                                        tablenames = c("bayes-pr2", "idtax-pr2"), 
                                        ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"))
  if (i == 1) {
    plot.list.pr2 <- list(r2way.pr2[[4]])
    plot.list2.pr2 <- list(a2way.pr2[[3]])
  }
  else {
    plot.list.pr2 <- c(plot.list.pr2, list(r2way.pr2[[4]]))
    plot.list2.pr2 <- c(plot.list2.pr2, list(a2way.pr2[[3]]))
  }
  
  x3 <- bayes.silva.list[[i]]
  x3[bayes.silva.conf < bootvec[i]] <- NA
  bayes.silva.list[[i]] <- x3
  
  x4 <- idtax.silva.list[[i]]
  x4[idtax.silva.conf < bootvec[i]] <- NA
  idtax.silva.list[[i]] <- x4
  
  # 2-way comparison of pr2 tables
  r2way.silva <- compare_byRank_2way(x3, x4,
                                   pltfilez = "none",
                                   tablenames = c("bayes-silva","idtax-silva"), 
                                   ranknamez = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))
  a2way.silva <- compare_assignments_2way(x3, x4,
                                     pltfilez = "none",
                                     tablenames = c("bayes-silva","idtax-silva"), 
                                     ranknamez = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))
  if (i == 1) {
    plot.list.silva <- list(r2way.silva[[4]])
    plot.list2.silva <- list(a2way.silva[[3]])
  }
  else {
    plot.list.silva <- c(plot.list.silva, list(r2way.silva[[4]]))
    plot.list2.silva <- c(plot.list2.silva, list(a2way.silva[[3]]))
  }
}

# # create table names and compare the resolution of all flavors of the pr2 taxonomy tables:
# bnam.pr2 <- sapply(bootvec.str, paste0, " bayes-pr2") # bayes pr2 table names...
# inam.pr2 <- sapply(bootvec.str, paste0, " idtax-pr2") # idtax pr2 names...
# tblnam <- c(bnam.pr2, inam.pr2)

# pr2 - bayes vs. idtax
rezcomp.bayes.pr2 <- compare_taxrez(bayes.pr2.list[[1]], bayes.pr2.list[[2]], bayes.pr2.list[[3]], bayes.pr2.list[[4]], bayes.pr2.list[[5]], 
                          pltfilez = "none",
                          tablenames = bootvec.str,
                          ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"))
rezcomp.idtax.pr2 <- compare_taxrez(idtax.pr2.list[[1]], idtax.pr2.list[[2]], idtax.pr2.list[[3]], idtax.pr2.list[[4]], idtax.pr2.list[[5]],
                                    pltfilez = "none",
                                    tablenames = bootvec.str,
                                    ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"))
# silva - bayes vs. idtax
rezcomp.bayes.silva <- compare_taxrez(bayes.silva.list[[1]], bayes.silva.list[[2]], bayes.silva.list[[3]], bayes.silva.list[[4]], bayes.silva.list[[5]], 
                                    pltfilez = "none",
                                    tablenames = bootvec.str,
                                    ranknamez = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))
rezcomp.idtax.silva <- compare_taxrez(idtax.silva.list[[1]], idtax.silva.list[[2]], idtax.silva.list[[3]], idtax.silva.list[[4]], idtax.silva.list[[5]],
                                    pltfilez = "none",
                                    tablenames = bootvec.str,
                                    ranknamez = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))

# stitch your plots together (with cowplot) and save
library("cowplot")
yl <- c(0,1) # y-limits
# pr2 comparisons by rank:
legend_b <- get_legend(plot.list.pr2[[1]])
p.2way.byRank.bootgradient.pr2 <- plot_grid(
  plot.list.pr2[[1]] + ggtitle(bootvec.str[1]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0,0.75,0.75,0.75), "cm")),
  plot.list.pr2[[2]] + ggtitle(bootvec.str[2]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0,0.75,0.75,0.75), "cm")),
  plot.list.pr2[[3]] + ggtitle(bootvec.str[3]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0,0.75,0.75,0.75), "cm")),
  plot.list.pr2[[4]] + ggtitle(bootvec.str[4]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  plot.list.pr2[[5]] + ggtitle(bootvec.str[5]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  legend_b,
  align = 'hv',
  labels = c("A.", "B.", "C.","D.", "E.",""),
  axis = 'l',
  hjust=-1,
  nrow=2
)
# adding a title:
title <- ggdraw() + 
  draw_label("bayes-pr2 vs. idtax-pr2: comparisons by rank",
    fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
p.2way.byRank.bootgradient.pr2 <- plot_grid(
  title, p.2way.byRank.bootgradient.pr2,
  ncol = 1,
  rel_heights = c(0.1, 1))
ggsave("bayesVidtax_boot_threshold/pr2_2way_byRank_boot40to80.pdf", p.2way.byRank.bootgradient.pr2, width = 15, height = 12, units = "in", device="pdf")

# pr2 comparisons by assignment:
p.2way.bootgradient.pr2 <- plot_grid(
  plot.list2.pr2[[1]] + ggtitle(bootvec.str[1]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  plot.list2.pr2[[2]] + ggtitle(bootvec.str[2]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  plot.list2.pr2[[3]] + ggtitle(bootvec.str[3]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  plot.list2.pr2[[4]] + ggtitle(bootvec.str[4]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  plot.list2.pr2[[5]] + ggtitle(bootvec.str[5]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  align = 'hv',
  labels = c("A.", "B.", "C.","D.", "E."),
  axis = 'l',
  hjust=-1,
  nrow=2
)
# adding a title:
title <- ggdraw() + 
  draw_label("bayes-pr2 vs. idtax-pr2: comparisons by assignment",
             fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
p.2way.bootgradient.pr2 <- plot_grid(
  title, p.2way.bootgradient.pr2,
  ncol = 1,
  rel_heights = c(0.1, 1))
ggsave("bayesVidtax_boot_threshold/pr2_2way_byAss_boot40to80.pdf", p.2way.bootgradient.pr2, width = 15, height = 12, units = "in", device="pdf")

# silva comparisons by rank:
legend_b <- get_legend(plot.list.silva[[1]])
p.2way.byRank.bootgradient.silva <- plot_grid(
  plot.list.silva[[1]] + ggtitle(bootvec.str[1]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  plot.list.silva[[2]] + ggtitle(bootvec.str[2]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  plot.list.silva[[3]] + ggtitle(bootvec.str[3]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  plot.list.silva[[4]] + ggtitle(bootvec.str[4]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  plot.list.silva[[5]] + ggtitle(bootvec.str[5]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  legend_b,
  align = 'hv',
  labels = c("A.", "B.", "C.","D.", "E.",""),
  axis = 'l',
  hjust=-1,
  nrow=2
)
# adding a title:
title <- ggdraw() + 
  draw_label("bayes-silva vs. idtax-silva: comparisons by rank",
             fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
p.2way.byRank.bootgradient.silva <- plot_grid(
  title, p.2way.byRank.bootgradient.silva,
  ncol = 1,
  rel_heights = c(0.1, 1))
ggsave("bayesVidtax_boot_threshold/silva_2way_byRank_boot40to80.pdf", p.2way.byRank.bootgradient.silva, width = 15, height = 12, units = "in", device="pdf")

# silva comparisons by assignment:
p.2way.bootgradient.silva <- plot_grid(
  plot.list2.silva[[1]] + ggtitle(bootvec.str[1]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  plot.list2.silva[[2]] + ggtitle(bootvec.str[2]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  plot.list2.silva[[3]] + ggtitle(bootvec.str[3]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  plot.list2.silva[[4]] + ggtitle(bootvec.str[4]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  plot.list2.silva[[5]] + ggtitle(bootvec.str[5]) + coord_cartesian(ylim = yl) + theme(legend.position="none",axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  align = 'hv',
  labels = c("A.", "B.", "C.","D.", "E."),
  axis = 'l',
  hjust=-1,
  nrow=2
)
# adding a title:
title <- ggdraw() + 
  draw_label("bayes-silva vs. idtax-silva: comparisons by assignment",
             fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
p.2way.bootgradient.silva <- plot_grid(
  title, p.2way.bootgradient.silva,
  ncol = 1,
  rel_heights = c(0.1, 1))
ggsave("bayesVidtax_boot_threshold/silva_2way_byAss_boot40to80.pdf", p.2way.bootgradient.silva, width = 15, height = 12, units = "in", device="pdf")

# resolution comparisons
p.rezcomps <- plot_grid(
  rezcomp.bayes.pr2[[2]] + ggtitle("bayes-pr2") + coord_cartesian(ylim = yl) + theme(axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  rezcomp.idtax.pr2[[2]] + ggtitle("idtax-pr2") + coord_cartesian(ylim = yl) + theme(axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  rezcomp.bayes.silva[[2]] + ggtitle("bayes-silva") + coord_cartesian(ylim = yl) + theme(axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  rezcomp.idtax.silva[[2]] + ggtitle("idtax-silva") + coord_cartesian(ylim = yl) + theme(axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
  align = 'hv',
  labels = c("A.", "B.", "C.","D."),
  axis = 'l',
  hjust=-1,
  nrow=2
)
ggsave("bayesVidtax_boot_threshold/all4_rezcomps_boot40to80.pdf", p.rezcomps, width = 10, height = 8, units = "in", device="pdf")

# preliminary pairwise comparisons of lca tax-tables:
rnam <- as.character(seq(from = 1, to = ncol(lca.pr2)-2, by = 1))
rnam <- sapply(rnam, FUN = function(x) {paste0("r",x)})
tblnam <- c("lca-pr2", "lca-silva")
# lca.comp.rank <- compare_byRank_2way(lca.pr2[,3:ncol(lca.pr2)], lca.silva[,3:ncol(lca.silva)],
#                                      pltfilez = "none",
#                                      tablenames = tblnam, 
#                                      ranknamez = rnam)
lca.comp.ass <- compare_assignments_2way(lca.pr2[,3:ncol(lca.pr2)], lca.silva[,3:ncol(lca.silva)],
                                         pltfilez = "none",
                                         tablenames = tblnam, 
                                         ranknamez = rnam)
p.lcacompz <- lca.comp.ass[[3]] + coord_cartesian(ylim = yl) + theme(axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm"))
# p.lcacompz <- plot_grid(
#   lca.comp.ass[[3]] + coord_cartesian(ylim = yl) + theme(axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
#   lca.comp.rank[[4]] + coord_cartesian(ylim = yl) + theme(axis.text.x = element_text(angle=45, hjust=1, size=12), plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm")),
#   align = 'hv',
#   labels = c("A.", "B.", "C.","D."),
#   axis = 'l',
#   hjust=-1,
#   nrow=1
# )
ggsave("bayesVidtax_boot_threshold/lca_unmapped_pr2vsSilva.pdf", p.lcacompz, width = 6, height = 4, units = "in", device="pdf")

# clear out the shit
rm(list=setdiff(ls(), c("bayes.pr2","bayes.pr2.conf","bayes.silva","bayes.silva.conf",
                        "idtax.pr2","idtax.pr2.conf","idtax.silva","idtax.silva.conf",
                        "lca.pr2", "lca.silva", "bootvec", "bootvec.str")))
# re-source your fcns:
source("~/Documents/R/amplicon_bioinformatics/package_deal/all_of_it.R")

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


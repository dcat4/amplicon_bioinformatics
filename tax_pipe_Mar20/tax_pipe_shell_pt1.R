# UNDER CONSTRUCTION

# up next:
# 1. I think the confidence dataframes might get shuffled when you merge since the merging sorts
# ... need to fix that... in bayes + idtax2df fcns...
# 2. add idtax2df_silva to helper_fcns so the changes propagate in both places

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

#### Step 1: Preliminary mapping of individual tax tables onto trait database 

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
lca.silva <- readRDS("initial_tax_tabs/LCA_pr2_rawdf_Mar20.rds")

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

identical(bayes.pr2$ASV, bayes.pr2$ASV)
identical(bayes.pr2$ASV, idtax.pr2$ASV)
identical(idtax.pr2$ASV, idtax.silva$ASV)  
identical(idtax.silva$ASV, lca.pr2$ASV)
identical(lca.pr2$ASV, lca.silva$ASV)

# let's peak at everything to make sure it's all good:
head(bayes.pr2)
head(bayes.silva)
head(idtax.pr2)
head(idtax.silva)
head(lca.pr2)
head(lca.silva)
# Noice.

# ok they look good. Remove ASV and svN col's from temp arrays of each
# little bootstrap threshold optimization:
bootvec <- seq(from = 10, to = 90, by = 10)
for (i in 1:length(bootvec)) {
  
}

# make sure boot-strapping thresholds are honored here...
# remove non-protists...

# then write a shell of my trait mapping and analysis functions

# preliminary pairwise comparisons of tax-tables too...

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


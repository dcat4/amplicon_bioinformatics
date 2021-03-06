# this is a script that takes the mock v9 data set and runs it through an "ensemble taxonomy workflow"
# compares the 3 individual tax tables to the expected, creates a couple different ensembles and compares those

# Kevin started this and I'm picking up where he left off -- I think I got it with a nice plot at the end, 
# but not the result I wanted/expected....

rm(list=ls())
setwd("~/Documents/R/amplicon_bioinformatics/mock_analysis/")
source("~/Documents/R/amplicon_bioinformatics/package_deal/all_of_it.R")

library("DECIPHER")
# rubric for aligning ASV numbers and sequences across datasets
rubber <- readDNAStringSet("mock_data/v9_mock_asvs_bothPrimers.fasta")

# read in and wrangle the individual tax tables:
lca.raw <- read.csv("mock_data/raw_LCA_out_mock_v9_pr2_May20.txt", header = FALSE, stringsAsFactors = FALSE)
lca <- LCA2df(lca.raw, rubber)

idtax.raw <- readRDS("mock_data/v9_mock_bothPrimers_idtax_pr2_0boot.rds")
idtax <- idtax2df_pr2(idtax.raw, boot=50, rubric=rubber)

bayes.raw <- readRDS("mock_data/v9_mock_bothPrimers_bayes_pr2_0boot.rds")
bayes <- bayestax2df(bayes.raw, boot=60, rubric=rubber)

# standardized table names for the pr2
table.names <- c("svN", "ASV", "kingdom", "supergroup", "division", "class", "order", "family", "genus", "species")
colnames(bayes) <- table.names
colnames(idtax) <- table.names

# setting up the parameters for mapping LCA
synonym.filepath <- "~/Documents/R/amplicon_bioinformatics/taxonomy_pipeline/tax_table_mapping/tax_synonyms_FINAL.csv"
# Bacteria and Archaea doesn't exist in pr2
nonexistent <- c('Bacteria', 'Archaea')
pr2 <- read.csv("~/Documents/R/amplicon_bioinformatics/taxonomy_pipeline/tax_table_mapping/pr2_all_tax.csv")
pr2 <- pr2[,-1]
lca.mapped <- taxmapper(taxin=lca, tax2map2=pr2,
                             exceptions=nonexistent,
                             synonym.file=synonym.filepath,
                             ignore.format=TRUE)
lca.mapped <- lca.mapped[[3]]
# sort the mapped output by svN
lca.mapped <- lca.mapped[order(nchar(as.character(lca.mapped$svN)),lca.mapped$svN), ]

exp.mock <- readRDS('mock_data/v9_mock_exptax_bothPrimers.rds')
colnames(exp.mock) <- table.names

make.NA <- function(x) if(is.character(x)||is.factor(x)){
  is.na(x) <- x=="NaN"; x} else {
    x}

exp.mock[] <- lapply(exp.mock, make.NA)

all.equal(nrow(bayes), nrow(idtax), nrow(lca.mapped))
all.equal(ncol(bayes), ncol(idtax), ncol(lca.mapped))

compare_results <- function(exp, ref, compare.cols) {
  # create a new data frame with the svN and ASV with the results
  results <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c('svN', 'ASV', 'exact', 'mis', 'over', 'under'))))
  
  # iterate through each row 
  # iterate from most general column to specific 
  nrows <- nrow(exp)
  for (row in 1:nrows) {
    n.same <- 0
    n.mis <- 0 
    n.over <- 0
    n.under <- 0
    for (col in compare.cols) {
      exp.tax <- exp[row, col]
      ref.tax <- ref[row, col]
      # over classified
      if (is.na(exp.tax) && !is.na(ref.tax)) {
        n.over <- n.over + 1
      }
      # under classified
      else if (!is.na(exp.tax) && is.na(ref.tax)) {
        n.under <- n.under + 1
      }
      # exact match
      else if ((is.na(exp.tax) && is.na(ref.tax)) || (exp.tax == ref.tax)) {
        n.same <- n.same + 1
      }
      # mis classified
      else {
        n.mis <- n.mis + 1
      }
    }
    r.row <- data.frame(matrix(rep(NA, 6), ncol=6, nrow = 1, dimnames=list(NULL, names(results))))
    r.row[, 'svN'] <- ref[row, 'svN']
    r.row[, 'ASV'] <- ref[row, 'ASV']
    r.row[, 'exact'] <- n.same
    r.row[, 'mis'] <- n.mis
    r.row[, 'over'] <- n.over
    r.row[, 'under'] <- n.under
    results <- rbind(results, r.row)
  }
  return(results)
}

bayes.mock.exp.c1 <- compare_results(exp.mock, bayes, table.names[-c(1:2)])
idtax.mock.exp.c1 <- compare_results(exp.mock, idtax, table.names[-c(1:2)])
lca.mock.exp.c1 <- compare_results(exp.mock, lca.mapped, table.names[-c(1:2)])

plot_results <- function(exp, ref, compare.cols) {
  n.same <- 0
  n.mis <- 0 
  n.over <- 0
  n.under <- 0
  
  # iterate through each row 
  # iterate from most general column to specific 
  nrows <- nrow(exp)
  for (row in 1:nrows) {
    curr.col<- 0
    for (col in compare.cols) {
      exp.tax <- exp[row, col]
      ref.tax <- ref[row, col]
      # over classified
      if (is.na(exp.tax) && !is.na(ref.tax)) {
        n.over <- n.over + 1
        break
      }
      # under classified
      else if (!is.na(exp.tax) && is.na(ref.tax)) {
        n.under <- n.under + 1
        break
      }
      # mis classified
      else if (!is.na(exp.tax) && !is.na(ref.tax) && (exp.tax != ref.tax)) {
        n.mis <- n.mis + 1
        break
      }
      # same
      else {
        curr.col <- curr.col + 1
      }
    }
    if (curr.col == length(compare.cols)) {
      n.same <- n.same + 1
    }
  }
  return(c(n.same, n.mis, n.over, n.under))
}

bayes.mock.exp.c2 <- plot_results(exp.mock, bayes, table.names[-c(1:2)])
idtax.mock.exp.c2 <- plot_results(exp.mock, idtax, table.names[-c(1:2)])
lca.mock.exp.c2 <- plot_results(exp.mock, lca.mapped, table.names[-c(1:2)])

# compute number of exact matches
nrow(bayes.mock.exp.c1[which(bayes.mock.exp.c1$exact == 8), ]) == bayes.mock.exp.c2[1]
nrow(idtax.mock.exp.c1[which(idtax.mock.exp.c1$exact == 8), ]) == idtax.mock.exp.c2[1]
nrow(lca.mock.exp.c1[which(lca.mock.exp.c1$exact == 8), ]) == lca.mock.exp.c2[1]

# compute number of mismatches
nrow(bayes.mock.exp.c1[which(bayes.mock.exp.c1$mis > 0), ]) == bayes.mock.exp.c2[2]
nrow(idtax.mock.exp.c1[which(idtax.mock.exp.c1$mis > 0), ]) == idtax.mock.exp.c2[2]
nrow(lca.mock.exp.c1[which(lca.mock.exp.c1$mis > 0), ]) == lca.mock.exp.c2[2]

# compute number of over
nrow(bayes.mock.exp.c1[which(bayes.mock.exp.c1$over > 0 & bayes.mock.exp.c1$mis == 0), ]) == bayes.mock.exp.c2[3]
nrow(idtax.mock.exp.c1[which(idtax.mock.exp.c1$over > 0 & idtax.mock.exp.c1$mis == 0), ]) == idtax.mock.exp.c2[3]
nrow(lca.mock.exp.c1[which(lca.mock.exp.c1$over > 0 & lca.mock.exp.c1$mis == 0), ]) == lca.mock.exp.c2[3]

# compute number of under
nrow(bayes.mock.exp.c1[which(bayes.mock.exp.c1$under > 0 & bayes.mock.exp.c1$mis == 0), ]) == bayes.mock.exp.c2[4]
nrow(idtax.mock.exp.c1[which(idtax.mock.exp.c1$under > 0 & idtax.mock.exp.c1$mis == 0), ]) == idtax.mock.exp.c2[4]
nrow(lca.mock.exp.c1[which(lca.mock.exp.c1$under > 0 & lca.mock.exp.c1$mis == 0), ]) == lca.mock.exp.c2[4]

# you need to add a warning to the consensus funciton if output is non-operational (like Na in the middle of an assignment, etc)
# compute the ensemble from the 3 tables
tblnam <- c("bayes-pr2", "idtax-pr2", "lca-pr2")
all.c <- consensus_tax_mostCom(bayes, idtax, lca.mapped, 
                               tablenames=tblnam, ranknamez=table.names,
                               tiebreakz=list(c("idtax-pr2", NA)), count.na=TRUE, weights=c(1,1,1))

# all.c.no.NA <- consensus_tax_mostCom(bayes, idtax, lca.mapped, 
#                                      tablenames=tblnam, ranknamez=table.names,
#                                      tiebreakz=list(c("idtax-pr2", NA)), count.na=FALSE, weights=c(1,1,1))

cbi <- consensus_tax_mostCom(bayes, idtax, 
                             tablenames=tblnam[1:2], ranknamez=table.names,
                              tiebreakz=list(c("idtax-pr2", NA)), count.na=TRUE, weights=c(1,1,1))
cbi.no.NA <- consensus_tax_mostCom(bayes, idtax,
                             tablenames=tblnam[1:2], ranknamez=table.names,
                             tiebreakz=list(c("idtax-pr2", NA)), count.na=FALSE, weights=c(1,1,1))

cbl <- consensus_tax_mostCom(bayes, lca.mapped, 
                             tablenames=tblnam[c(1,3)], ranknamez=table.names,
                             tiebreakz=list(c("bayes-pr2", NA)), count.na=TRUE, weights=c(1,1,1))
cbl.no.NA <- consensus_tax_mostCom(bayes, lca.mapped,
                                   tablenames=tblnam[c(1,3)], ranknamez=table.names,
                                   tiebreakz=list(c("bayes-pr2", NA)), count.na=FALSE, weights=c(1,1,1))
cil <- consensus_tax_mostCom(idtax, lca.mapped, 
                             tablenames=tblnam[2:3], ranknamez=table.names,
                             tiebreakz=list(c("idtax-pr2", NA)), count.na=TRUE, weights=c(1,1,1))
cil.no.NA <- consensus_tax_mostCom(idtax, lca.mapped, 
                                   tablenames=tblnam[2:3], ranknamez=table.names,
                                   tiebreakz=list(c("idtax-pr2", NA)), count.na=FALSE, weights=c(1,1,1))

cbi.compare <- plot_results(exp.mock, cbi, table.names[-c(1:2)])
cbi.no.NA.compare <- plot_results(exp.mock, cbi.no.NA, table.names[-c(1:2)])
cbl.compare <- plot_results(exp.mock, cbl, table.names[-c(1:2)])
cbl.no.NA.compare <- plot_results(exp.mock, cbl.no.NA, table.names[-c(1:2)])
cil.compare <- plot_results(exp.mock, cil, table.names[-c(1:2)])
cil.no.NA.compare <- plot_results(exp.mock, cil.no.NA, table.names[-c(1:2)])
all.c.compare <- plot_results(exp.mock, all.c, table.names[-c(1:2)])

# create a data frame
vals <- rbind(bayes.mock.exp.c2, idtax.mock.exp.c2, lca.mock.exp.c2, 
              cbi.compare, cbi.no.NA.compare, cbl.compare, cbl.no.NA.compare, 
              cil.compare, cil.no.NA.compare, all.c.compare)
vals

dn <- list(c('bayes', 'idtax', 'lca', 'bayes+idtax', 'bayes+idtax-noNA', 'bayes+lca', 'bayes+lca-noNA', 'idtax+lca', 'idtax+lca-noNA', 'all 3'), 
           c('exact', 'mis', 'over', 'under'))
df <- data.frame(matrix(vals, ncol=ncol(vals), nrow = nrow(vals), dimnames=dn))
df <- (df / rowSums(df)) * 100

bloop <- dn[[1]]
table <- c()
for (i in 1:length(bloop)) {
  table <- append(table, rep(bloop[i], 4))
}
class <- c(rep(c('exact', 'mis', 'over', 'under'), length(bloop)))

library("reshape2")
df$table <- rownames(df)
data <- melt(df)

library("ggplot2")
ggplot(data, aes(fill=variable, y=value, x=table)) +
  geom_bar(position=position_dodge(width=0.8), stat='identity') +
  geom_text(aes(label=round(value, digits = 2)), position=position_dodge(width=0.8), vjust=-0.25, size=6) +
  labs(x="taxonomy table", y='% of ASVs') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 25), axis.title.x = element_text(size = 25, face="bold"),
        axis.text.y = element_text(size = 25), axis.title.y = element_text(size = 25, face="bold"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white"),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black")) + 
  scale_fill_discrete(name = "classification relative to expected") +
  theme(legend.text = element_text(size=20),
        legend.title = element_text(size=20)  )

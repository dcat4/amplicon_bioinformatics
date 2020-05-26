# this is a shell script that will create ensemble taxonomies and compare individual taxtabs
# to try to start quantifying uncertainties

# got it set up and works pretty nice :)

rm(list=ls())
setwd("~/Documents/R/amplicon_bioinformatics/tax_pipe_Mar20/")
source("~/Documents/R/amplicon_bioinformatics/package_deal/all_of_it.R")
xx <- readRDS(file = "all_mapped_taxtabs_protistOnly.rds")

seqtab <- readRDS("seqtab_nochime_Mar20.rds")
samsubs <- readRDS("sample_subsets.rds")

# construct consensus tax based on best resolution, w/ LCA-like as a tie-breaker
ctax.br <- consensus_tax_bestRez(xx[[1]], xx[[2]],xx[[3]], xx[[4]],xx[[5]], xx[[6]],
                      tablenames = names(xx), 
                      ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"),
                      tiebreakz = "LCAlike")

ctax.br <- ctax.br[[1]] # extract the consensus array
# this bit copy/pasted, just need to tweak for what I want here...
for (i in 1:length(xx)) {

  x1 <-xx[[i]]
  n1 <- names(xx)[i]

  r2way.br <- compare_byRank_2way(ctax.br, x1,
                               pltfilez = "none",
                               tablenames = c("en-best-rez", n1), 
                               ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"))
  a2way.br <- compare_assignments_2way(ctax.br, x1,
                                    pltfilez = "none",
                                    tablenames = c("en-best-rez", n1), 
                                    ranknamez = c("Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"))
  if (i == 1) {
    plot.list.br <- list(r2way.br[[4]])
    plot.list2.br <- list(a2way.br[[3]])
  }
  else {
    plot.list.br <- c(plot.list.br, list(r2way.br[[4]]))
    plot.list2.br <- c(plot.list2.br, list(a2way.br[[3]]))
  }
}


# stitch the plots together w/ cowplot
library("cowplot")
yl <- c(0,1) # y-limits

p.byrank <- plot_grid(
  plot.list.br[[1]] + coord_cartesian(ylim = yl) + ggtitle(names(xx)[1]),
  plot.list.br[[2]] + coord_cartesian(ylim = yl) + ggtitle(names(xx)[2]),
  plot.list.br[[3]] + coord_cartesian(ylim = yl) + ggtitle(names(xx)[3]),
  plot.list.br[[4]] + coord_cartesian(ylim = yl) + ggtitle(names(xx)[4]),
  plot.list.br[[5]] + coord_cartesian(ylim = yl) + ggtitle(names(xx)[5]),
  plot.list.br[[6]] + coord_cartesian(ylim = yl) + ggtitle(names(xx)[6]),
  align = 'hv',
  labels = c("A.", "B.", "C.","D.","E.","F."),
  axis = 'l',
  hjust=-1,
  nrow = 3
)
ggsave("comp_results_taxtab2ensembles/all6tab_vs_ensembleBestRez.pdf", p.byrank, width = 15, height = 20, units = "in", device="pdf")

p.byass <- plot_grid(
  plot.list2.br[[1]] + coord_cartesian(ylim = yl) + ggtitle(names(xx)[1]),
  plot.list2.br[[2]] + coord_cartesian(ylim = yl) + ggtitle(names(xx)[2]),
  plot.list2.br[[3]] + coord_cartesian(ylim = yl) + ggtitle(names(xx)[3]),
  plot.list2.br[[4]] + coord_cartesian(ylim = yl) + ggtitle(names(xx)[4]),
  plot.list2.br[[5]] + coord_cartesian(ylim = yl) + ggtitle(names(xx)[5]),
  plot.list2.br[[6]] + coord_cartesian(ylim = yl) + ggtitle(names(xx)[6]),
  align = 'hv',
  labels = c("A.", "B.", "C.","D.","E.","F."),
  axis = 'l',
  hjust=-1,
  nrow = 3
)
ggsave("comp_results_taxtab2ensembles/all6tab_vs_ensembleBestRez2.pdf", p.byass, width = 15, height = 20, units = "in", device="pdf")




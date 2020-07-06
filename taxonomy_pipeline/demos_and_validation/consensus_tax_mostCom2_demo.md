consensus\_tax\_mostCom2 Demo
================
D Catlett
7/2/2020

Overview
--------

This is a demo showing the functionality of the re-written *consensus\_tax\_mostCom2* algorithm that creates consensus taxonomic assignments from multiple input taxonomy tables.

### consensus\_tax\_mostCom Algorithm Description

consensus\_tax\_mostCom generates a consensus taxonomic assignment for each ASV/OTU in a data set by using the taxonomic assignment found in a majority (&gt;50%) of independent input taxonomic assignment data sets (produced by different algorithms, reference databases, or both). The algorithm iterates from the root to the species rank and assigns taxonomic names found in a majority of input taxonomy tables. Non-assignments indicated by NA are considered in this calculation by default.

The algorithm offers several user-controllable parameters that make it's behavior more flexible. These are demonstrated below and include options to:

1.  Relax the majority threshold to a user-specified value such that consensus assignments are determined by finding the taxonomic name with a frequency higher than a user-specified threshold.

2.  Prioritize assignments of particular input taxonomy tables in the event where multiple taxonomic names have equivalent, maximal frequencies across the input taxonomy tables.

3.  Weight the assignments of particular input taxonomy tables. Weights are applied in the frequency calculations. This is independent of the tie-breaking option in \#2.

4.  Ignore non-assignemnts (NA) in the consensus calculation.

### Setting Things Up:

To get started, we'll clear out our environment, load in the reshape2 package ofr later, set our wd, and read in taxonomy tables: The taxonomy tables used here come from implementations of the RDP Bayesian classifer, the new idtaxa algorithm, and MEGAN'S LCA algorithm against the pr2 reference databases. Our amplicon data set is an 18S-V9 tag sequencing project from the coastal ocean.

You can do this with any set of taxonomy tables assuming you format them properly. To follow along with this demo, grab the taxonomy tables in the directory listed below and run all the code you see below.

We'll use the R packages dplyr and DECIPHER to get through this data, so go ahead and install those too.

``` r
rm(list = ls())
library("dplyr")
library("DECIPHER")

# setwd and read in your datasets:
setwd("~/Documents/R/amplicon_bioinformatics/tax_pipe_Mar20/")
dd <- "~/Documents/R/amplicon_bioinformatics/package_deal/"

# read in all our functions:
source(paste0(dd,"all_of_it.R"))

# read in  the data
bayes.pr2 <- readRDS("initial_tax_tabs/bayes_pr2_0boot_Mar20.rds")
idtax.pr2 <- readRDS("initial_tax_tabs/idtax_pr2_0boot_Mar20.rds")
lca.pr2 <- readRDS("initial_tax_tabs/LCA_pr2_rawdf_Mar20.rds")
```

### Arranging and formating our taxonomy tables for running the algorithm:

Here we'll reformat our taxonomy tables using some homebrewed helper functions, map the lca table (with NCBI's taxonomic nomenclature) onto pr2's taxonomic nomenclature using our *taxmapper* functions, and sort each taxonomy table by ASV sequences so they're all in the same order.

``` r
rubber <- readDNAStringSet(filepath = "~/Documents/R/amplicon_bioinformatics/tax_pipe_Mar20/blaster/allASVs4blast.fasta")
# convert idtax and bayes to dataframes - keep boot = 0  and return conf for now...
# note that for idtaxa, order of ASVs must be retrieved from dada2 seqtab - rubric represents that order
# and adds ASV names and sequences to output dataframe:
xx <- bayestax2df(bayes.pr2, boot = 60, rubric = rubber, return.conf = TRUE)
bayes.pr2 <- xx[[1]]

xx <- idtax2df_pr2(idtax.pr2, boot = 50, rubric = rubber, return.conf = TRUE)
idtax.pr2 <- xx[[1]]

# re-order your 6 tables by sorting according to ASV:
ii <- base::sort(bayes.pr2$ASV, index.return = TRUE)
bayes.pr2 <- bayes.pr2[ii$ix,]
ii <- base::sort(idtax.pr2$ASV, index.return = TRUE)
idtax.pr2 <- idtax.pr2[ii$ix,]

tt <- read.csv(file = "~/Documents/R/amplicon_bioinformatics/taxonomy_pipeline/tax_table_mapping/pr2_all_tax.csv", stringsAsFactors = FALSE)
tt <- tt[, -c(1)]
# map lca.pr2 to pr2 taxonomy:
lca.pr2 <- taxmapper(lca.pr2, tax2map2 = tt, 
                     ignore.format = TRUE, exceptions = c("Bacteria", "Archaea"),
                     synonym.file = "~/Documents/R/amplicon_bioinformatics/taxonomy_pipeline/tax_table_mapping/tax_synonyms_FINAL.csv", outfilez = "none")
lca.pr2 <- lca.pr2[[3]]
ii <- base::sort(lca.pr2$ASV, index.return = TRUE)
lca.pr2 <- lca.pr2[ii$ix,]
```

You can run this for a sanity check:

``` r
# check that they're all in the same order
identical(bayes.pr2$ASV, idtax.pr2$ASV)
identical(bayes.pr2$ASV, lca.pr2$ASV)
```

...and this to see what the data sets look like. These data sets are available in the test-data directory.

``` r
head(bayes.pr2)
head(idtax.pr2)
head(lca.pr2)
```

The algorithm expects that the ASV's in each taxonomy table are all in the same order and employ a common taxonomic nomenclature (naming and ranking convention). If your tables don't meet these requirements, we recommend standardizing taxonomic nomenclatures using our taxonomy mapping algorithm *taxmapper* and sorting the rows of each of your taxonomy tables according to the ASV sequences in your data set (see above for demo's).

### First run

We've pre-selected a few ASVs to illustrate some test cases. Here we'll set up our data to make the function calls a little more streamlined, then step through a series of algorithm runs and check the results.

``` r
# set up the data:
x <- list(bayes.pr2, idtax.pr2, lca.pr2)
names(x) <- c("bayes.pr2", "idtax.pr2", "lca.pr2")
nn <- c("svN", "ASV", "kingdom", "supergroup", "division","class","order","family","genus","species")
for (i in 1:length(x)){
  bloop <- x[[i]]
  colnames(bloop) <- nn
  x[[i]] <- bloop
}
# subset test cases
testers <- c("sv20747", "sv20858", "sv14136", "sv17278", "sv20738","sv3579", "sv20673", "sv8382", "sv4298")
x <- lapply(x, function(y) y[y$svN %in% testers , ])
# run the algorithm
ctax <- consensus_tax_mostCom2(x[[1]], x[[2]], x[[3]], 
                               tablenames = names(x), ranknamez = c("kingdom", "supergroup", "division","class","order","family","genus","species"),
                               tiebreakz = "none", count.na=TRUE, trueMajority=TRUE, threshold = 0.5, weights=rep(1, length(x)))
```

    ## no Frankenstein assignments detected

``` r
# compile results for easy display later
test.results <- rbind(cbind(data.frame(tbl = c(rep(names(x)[1], nrow(ctax))), stringsAsFactors = FALSE), x[[1]]),
  cbind(data.frame(tbl = c(rep(names(x)[2], nrow(ctax))), stringsAsFactors = FALSE), x[[2]]),
  cbind(data.frame(tbl = c(rep(names(x)[3], nrow(ctax))), stringsAsFactors = FALSE), x[[3]]),
  cbind(data.frame(tbl = c(rep("consensus", nrow(ctax))), stringsAsFactors = FALSE), ctax))
# helper function for subsetting particular test cases in the results:
subsetter <- function(y, sv) {
  y[y$svN == sv , ]
}
# Look at "sv20747" first:
subsetter(test.results, "sv20747")
```

    ##             tbl     svN
    ## 3353  bayes.pr2 sv20747
    ## 20747 idtax.pr2 sv20747
    ## 13251   lca.pr2 sv20747
    ## 7     consensus sv20747
    ##                                                                                                                                 ASV
    ## 3353  GCACCCACCGATTGAAAAGCCCGGTGAAGAATCGGGATTGTAGCGTTGTCCTTCATTGGACATTGCCGTGAGAACCTTTCTGAACCTTGTTTTTTAGAGGAAGGTGAAGTCGTAACAAGGTCTCT
    ## 20747 GCACCCACCGATTGAAAAGCCCGGTGAAGAATCGGGATTGTAGCGTTGTCCTTCATTGGACATTGCCGTGAGAACCTTTCTGAACCTTGTTTTTTAGAGGAAGGTGAAGTCGTAACAAGGTCTCT
    ## 13251 GCACCCACCGATTGAAAAGCCCGGTGAAGAATCGGGATTGTAGCGTTGTCCTTCATTGGACATTGCCGTGAGAACCTTTCTGAACCTTGTTTTTTAGAGGAAGGTGAAGTCGTAACAAGGTCTCT
    ## 7     GCACCCACCGATTGAAAAGCCCGGTGAAGAATCGGGATTGTAGCGTTGTCCTTCATTGGACATTGCCGTGAGAACCTTTCTGAACCTTGTTTTTTAGAGGAAGGTGAAGTCGTAACAAGGTCTCT
    ##         kingdom    supergroup   division           class             order
    ## 3353  Eukaryota Stramenopiles Ochrophyta Bacillariophyta Bacillariophyta_X
    ## 20747 Eukaryota Stramenopiles Ochrophyta Bacillariophyta Bacillariophyta_X
    ## 13251 Eukaryota Stramenopiles Ochrophyta Bacillariophyta Bacillariophyta_X
    ## 7     Eukaryota Stramenopiles Ochrophyta Bacillariophyta Bacillariophyta_X
    ##                           family       genus species
    ## 3353  Polar-centric-Mediophyceae Chaetoceros    <NA>
    ## 20747 Polar-centric-Mediophyceae Chaetoceros    <NA>
    ## 13251 Polar-centric-Mediophyceae Chaetoceros    <NA>
    ## 7     Polar-centric-Mediophyceae Chaetoceros    <NA>

Now we have our first glimpse of the algorithms output. A message is displayed telling us no Frankenstein assignments were detected. This is just an internal QC and so this is a good thing.

We also see that for this test case, all 3 input taxonomies perfectly agreed in their taxonomic assignments, so the consensus assignment is identical to all 3. Looking at another:

``` r
subsetter(test.results, "sv3579")
```

    ##             tbl    svN
    ## 2867  bayes.pr2 sv3579
    ## 3579  idtax.pr2 sv3579
    ## 15624   lca.pr2 sv3579
    ## 5     consensus sv3579
    ##                                                                                                                                      ASV
    ## 2867  ACTCCTACCAATTGAATGATCCATGAAGTGTTTGGATTACATTGAAGATGGTGGTTTGCCGCTGTCGACGTCATGAGAAGTTCATTGAACCTTATCATTTAGAGGAAGGAGAAGTCATAACAAGGTTACC
    ## 3579  ACTCCTACCAATTGAATGATCCATGAAGTGTTTGGATTACATTGAAGATGGTGGTTTGCCGCTGTCGACGTCATGAGAAGTTCATTGAACCTTATCATTTAGAGGAAGGAGAAGTCATAACAAGGTTACC
    ## 15624 ACTCCTACCAATTGAATGATCCATGAAGTGTTTGGATTACATTGAAGATGGTGGTTTGCCGCTGTCGACGTCATGAGAAGTTCATTGAACCTTATCATTTAGAGGAAGGAGAAGTCATAACAAGGTTACC
    ## 5     ACTCCTACCAATTGAATGATCCATGAAGTGTTTGGATTACATTGAAGATGGTGGTTTGCCGCTGTCGACGTCATGAGAAGTTCATTGAACCTTATCATTTAGAGGAAGGAGAAGTCATAACAAGGTTACC
    ##         kingdom     supergroup     division         class           order
    ## 2867  Eukaryota Archaeplastida Streptophyta          <NA>            <NA>
    ## 3579       <NA>           <NA>         <NA>          <NA>            <NA>
    ## 15624 Eukaryota Archaeplastida Streptophyta Embryophyceae Embryophyceae_X
    ## 5     Eukaryota Archaeplastida Streptophyta          <NA>            <NA>
    ##                 family genus species
    ## 2867              <NA>  <NA>    <NA>
    ## 3579              <NA>  <NA>    <NA>
    ## 15624 Embryophyceae_XX Pinus    <NA>
    ## 5                 <NA>  <NA>    <NA>

In this scenario, because we set *count.na = TRUE*, we see that the consensus assignment is an unassigned Streptophyta. The lca.pr2 assignment was more resolved than this, but the overall majority assignment only extends to the Division rank. One more:

``` r
subsetter(test.results, "sv4298")
```

    ##             tbl    svN
    ## 4812  bayes.pr2 sv4298
    ## 4298  idtax.pr2 sv4298
    ## 13772   lca.pr2 sv4298
    ## 9     consensus sv4298
    ##                                                                                                                                ASV
    ## 4812  GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ## 4298  GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ## 13772 GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ## 9     GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ##         kingdom    supergroup   division            class              order
    ## 4812  Eukaryota Stramenopiles   Opalozoa           MAST-3            MAST-3B
    ## 4298  Eukaryota Stramenopiles       <NA>             <NA>               <NA>
    ## 13772 Eukaryota Stramenopiles Ochrophyta Dictyochophyceae Dictyochophyceae_X
    ## 9     Eukaryota Stramenopiles       <NA>             <NA>               <NA>
    ##             family      genus        species
    ## 4812     MAST-3B_X MAST-3B_XX MAST-3B_XX_sp.
    ## 4298          <NA>       <NA>           <NA>
    ## 13772 Pedinellales       <NA>           <NA>
    ## 9             <NA>       <NA>           <NA>

Here, the consensus assignment is an unidentified Stramenopile, which is really all we can infer from these input taxonomy tables since they each have a different assignment at the Division rank. Now we'll explore the results for this ASV as we alter certain parameters in the consensus taxonomy algorithm...

### Second run

Here we'll run the algorithm again and look at our *sv4298* above, but we'll set *trueMajority=FALSE* so we don't need a name to occur &gt;50% of the time for it to be assigned to the consensus. We'll alter the *threshold* required for the consensus assignment to *0.3* so that only one instance is needed for the name to be assigned. And finally, we'll specify a tie-breaker with *tiebreakz = list(c("bayes.pr2", NA))*. Altogether, this tells the algorithm to prioritize all assignments by the Bayes-pr2 taxonomy table in the event that those assignments occur at a proportion &gt; 0.3 across all taxonomy tables and are tied for the maximal frequency of the assignments. This is a little redundant in this example but you could imagine a scenario with more input taxonomy tables where this would make more sense.

``` r
# run the algorithm
ctax <- consensus_tax_mostCom2(x[[1]], x[[2]], x[[3]], 
                               tablenames = names(x), ranknamez = c("kingdom", "supergroup", "division","class","order","family","genus","species"),
                               tiebreakz = list(c("bayes.pr2", NA)), count.na=TRUE, trueMajority=FALSE, threshold = 0.3, weights=rep(1, length(x)))
```

    ## no Frankenstein assignments detected

``` r
# compile results for easy display later
test.results <- rbind(cbind(data.frame(tbl = c(rep(names(x)[1], nrow(ctax))), stringsAsFactors = FALSE), x[[1]]),
  cbind(data.frame(tbl = c(rep(names(x)[2], nrow(ctax))), stringsAsFactors = FALSE), x[[2]]),
  cbind(data.frame(tbl = c(rep(names(x)[3], nrow(ctax))), stringsAsFactors = FALSE), x[[3]]),
  cbind(data.frame(tbl = c(rep("consensus", nrow(ctax))), stringsAsFactors = FALSE), ctax))

subsetter(test.results, "sv4298")
```

    ##             tbl    svN
    ## 4812  bayes.pr2 sv4298
    ## 4298  idtax.pr2 sv4298
    ## 13772   lca.pr2 sv4298
    ## 9     consensus sv4298
    ##                                                                                                                                ASV
    ## 4812  GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ## 4298  GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ## 13772 GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ## 9     GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ##         kingdom    supergroup   division            class              order
    ## 4812  Eukaryota Stramenopiles   Opalozoa           MAST-3            MAST-3B
    ## 4298  Eukaryota Stramenopiles       <NA>             <NA>               <NA>
    ## 13772 Eukaryota Stramenopiles Ochrophyta Dictyochophyceae Dictyochophyceae_X
    ## 9     Eukaryota Stramenopiles   Opalozoa           MAST-3            MAST-3B
    ##             family      genus        species
    ## 4812     MAST-3B_X MAST-3B_XX MAST-3B_XX_sp.
    ## 4298          <NA>       <NA>           <NA>
    ## 13772 Pedinellales       <NA>           <NA>
    ## 9        MAST-3B_X       <NA>           <NA>

So now our consensus assignment matches the bayes input taxonomy table down to the family rank. This is because we had 3-way ties across the 3 input taxonomy tables at the division, class, and order ranks, and because we specified the bayes table as our tie-breaker, those assignments were used in the consensus. At the genus rank however, the non-assignments outnumber the bayes assignment 2:1. This means no tie-breaking was required, and the consensus is not assigned (or assigned NA).

### Third run

In this next run, we'll remove the tie-breaker, but we'll also set *count.na=FALSE*. This tells the algorithm not to consider NA's as an assignment. This would be used in situations where you can be a little less conservative in protecting against over-assignments. It's important to note that this does not necessarily result in no NA assignments in the consensus as there still must be a most-common assignment present to be included in the consensus.

``` r
# run the algorithm
ctax <- consensus_tax_mostCom2(x[[1]], x[[2]], x[[3]], 
                               tablenames = names(x), ranknamez = c("kingdom", "supergroup", "division","class","order","family","genus","species"),
                               tiebreakz = "none", count.na=FALSE, trueMajority=FALSE, threshold = 0.3, weights=rep(1, length(x)))
```

    ## no Frankenstein assignments detected

``` r
# compile results for easy display later
test.results <- rbind(cbind(data.frame(tbl = c(rep(names(x)[1], nrow(ctax))), stringsAsFactors = FALSE), x[[1]]),
  cbind(data.frame(tbl = c(rep(names(x)[2], nrow(ctax))), stringsAsFactors = FALSE), x[[2]]),
  cbind(data.frame(tbl = c(rep(names(x)[3], nrow(ctax))), stringsAsFactors = FALSE), x[[3]]),
  cbind(data.frame(tbl = c(rep("consensus", nrow(ctax))), stringsAsFactors = FALSE), ctax))

subsetter(test.results, "sv4298")
```

    ##             tbl    svN
    ## 4812  bayes.pr2 sv4298
    ## 4298  idtax.pr2 sv4298
    ## 13772   lca.pr2 sv4298
    ## 9     consensus sv4298
    ##                                                                                                                                ASV
    ## 4812  GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ## 4298  GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ## 13772 GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ## 9     GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ##         kingdom    supergroup   division            class              order
    ## 4812  Eukaryota Stramenopiles   Opalozoa           MAST-3            MAST-3B
    ## 4298  Eukaryota Stramenopiles       <NA>             <NA>               <NA>
    ## 13772 Eukaryota Stramenopiles Ochrophyta Dictyochophyceae Dictyochophyceae_X
    ## 9     Eukaryota Stramenopiles       <NA>             <NA>               <NA>
    ##             family      genus        species
    ## 4812     MAST-3B_X MAST-3B_XX MAST-3B_XX_sp.
    ## 4298          <NA>       <NA>           <NA>
    ## 13772 Pedinellales       <NA>           <NA>
    ## 9             <NA>       <NA>           <NA>

And we see that our consensus is only assigned to the supergroup rank because there was a 1-1 tie at all lower ranks between the bayes and lca tables. Because we did not specify a tie-breaker or weight one of these tables, there was no most common taxonomic assignment that could be assigned as the consensus.

### Fourth run

So as you probably guessed, here we're weighting the bayes table double that of the other 2, and still not counting NA's. Let's see what that does:

``` r
# run the algorithm
ctax <- consensus_tax_mostCom2(x[[1]], x[[2]], x[[3]], 
                               tablenames = names(x), ranknamez = c("kingdom", "supergroup", "division","class","order","family","genus","species"),
                               tiebreakz = "none", count.na=FALSE, trueMajority=FALSE, threshold = 0.3, weights=c(2, 1, 1))
```

    ## no Frankenstein assignments detected

``` r
# compile results for easy display later
test.results <- rbind(cbind(data.frame(tbl = c(rep(names(x)[1], nrow(ctax))), stringsAsFactors = FALSE), x[[1]]),
  cbind(data.frame(tbl = c(rep(names(x)[2], nrow(ctax))), stringsAsFactors = FALSE), x[[2]]),
  cbind(data.frame(tbl = c(rep(names(x)[3], nrow(ctax))), stringsAsFactors = FALSE), x[[3]]),
  cbind(data.frame(tbl = c(rep("consensus", nrow(ctax))), stringsAsFactors = FALSE), ctax))

subsetter(test.results, "sv4298")
```

    ##             tbl    svN
    ## 4812  bayes.pr2 sv4298
    ## 4298  idtax.pr2 sv4298
    ## 13772   lca.pr2 sv4298
    ## 9     consensus sv4298
    ##                                                                                                                                ASV
    ## 4812  GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ## 4298  GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ## 13772 GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ## 9     GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ##         kingdom    supergroup   division            class              order
    ## 4812  Eukaryota Stramenopiles   Opalozoa           MAST-3            MAST-3B
    ## 4298  Eukaryota Stramenopiles       <NA>             <NA>               <NA>
    ## 13772 Eukaryota Stramenopiles Ochrophyta Dictyochophyceae Dictyochophyceae_X
    ## 9     Eukaryota Stramenopiles   Opalozoa           MAST-3            MAST-3B
    ##             family      genus        species
    ## 4812     MAST-3B_X MAST-3B_XX MAST-3B_XX_sp.
    ## 4298          <NA>       <NA>           <NA>
    ## 13772 Pedinellales       <NA>           <NA>
    ## 9        MAST-3B_X MAST-3B_XX MAST-3B_XX_sp.

Now our consensus is assigned all the way to species! This is because at the ranks where idtax is NA, the bayes name is weighted 2x so the bayes assignment is counted twice and is thus the consensus, despite the fact that the LCA assignment is different.

Fifth run
---------

Now we'll remove the weights and prioritize the lca table just to give it some time to shine:

``` r
# run the algorithm
ctax <- consensus_tax_mostCom2(x[[1]], x[[2]], x[[3]], 
                               tablenames = names(x), ranknamez = c("kingdom", "supergroup", "division","class","order","family","genus","species"),
                               tiebreakz = list(c("lca.pr2", NA)), count.na=FALSE, trueMajority=FALSE, threshold = 0.3, weights=c(1, 1, 1))
```

    ## no Frankenstein assignments detected

``` r
# compile results for easy display later
test.results <- rbind(cbind(data.frame(tbl = c(rep(names(x)[1], nrow(ctax))), stringsAsFactors = FALSE), x[[1]]),
  cbind(data.frame(tbl = c(rep(names(x)[2], nrow(ctax))), stringsAsFactors = FALSE), x[[2]]),
  cbind(data.frame(tbl = c(rep(names(x)[3], nrow(ctax))), stringsAsFactors = FALSE), x[[3]]),
  cbind(data.frame(tbl = c(rep("consensus", nrow(ctax))), stringsAsFactors = FALSE), ctax))

subsetter(test.results, "sv4298")
```

    ##             tbl    svN
    ## 4812  bayes.pr2 sv4298
    ## 4298  idtax.pr2 sv4298
    ## 13772   lca.pr2 sv4298
    ## 9     consensus sv4298
    ##                                                                                                                                ASV
    ## 4812  GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ## 4298  GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ## 13772 GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ## 9     GCACCTACCGATTGAATGGTCCGGTGAGATCTTCGGACTGCAGCGAAAGTCAGCAATGAGTTAGTCGCGGAAAGTTGATCAAACCTTACCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC
    ##         kingdom    supergroup   division            class              order
    ## 4812  Eukaryota Stramenopiles   Opalozoa           MAST-3            MAST-3B
    ## 4298  Eukaryota Stramenopiles       <NA>             <NA>               <NA>
    ## 13772 Eukaryota Stramenopiles Ochrophyta Dictyochophyceae Dictyochophyceae_X
    ## 9     Eukaryota Stramenopiles Ochrophyta Dictyochophyceae Dictyochophyceae_X
    ##             family      genus        species
    ## 4812     MAST-3B_X MAST-3B_XX MAST-3B_XX_sp.
    ## 4298          <NA>       <NA>           <NA>
    ## 13772 Pedinellales       <NA>           <NA>
    ## 9     Pedinellales       <NA>           <NA>

And we see that the algorithm has correctly broken the ties. Sickarooni!!!

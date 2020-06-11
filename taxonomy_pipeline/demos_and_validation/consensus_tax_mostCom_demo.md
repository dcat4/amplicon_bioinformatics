consensus\_tax\_mostCom Demo
================
K Son
April 30, 2020

## Overview

Here, we step through implementations of the consensus\_tax\_mostCom.R
algorithm to demonstrate proper uses and anticipated results. This
algorithm is built to utilize information from multiple taxonomy tables
(defined here as a table of ASVs and correspnding taxonomic assignments)
in order to improve the resolution and/or accuracy of taxonomic
annotations of amplicon sequences. It incorporates information from
multiple taxonomy tables and determines a “consensus taxonomy” for each
ASV in your data set. The algorithm requires that all taxonomy tables
follow the same taxonomic naming and ranking conventions, that the order
of columns in each taxonomy table follows the taxonomic ranking
heirarchy, and that the order of rows (ASVs) in each of the input
taxonomy tables is the same. If these rules are not followed, spurious
results are likely.

### consensus\_tax\_mostCom Algorithm Description

consensus\_tax\_mostCom generates a consensus taxonomy by using the most
frequent taxonomy starting at the most specific ranking across all
user-specified input taxonomy tables. In other words, it’s designed to
get the the taxonomic annotation where it is the majority among other
taxonomic assignments. Unresolved taxonomic assignments in each taxonomy
table should be indicated by NA as there is a parameter for the
algorithm to consider NA as a taxonomic assignment or not. There will be
occasions when multiple input taxonomy tables are ‘tied’ for the most
frequeny taxonomy for a given ASV; thus, after assigning a consensus
taxonomy to the ASVs where a single taxonomy provides the most common
taxonomy, the algorithm uses a series of user-speicifed rules to break
the remaining ties. Rules you may specify are demonstarted below.

### Start ’er up:

We’ll clear out our environment, load in the reshape2 package ofr later,
set our wd, and read in taxonomy tables: The taxonomy tables used here
come from implementations of the RDP Bayesian classifer, the new idtaxa
algorithm, and MEGAN’S LCA algorithm against both the Silva and pr2
reference databases. Our amplicon data set is an 18S-V9 tag sequencing
project from the coastal ocean.

You can do this with any taxonomy tables assuming you format them
properly. To follow along with this demo, grab the taxonomy tables in
the “test\_data” directory of this repository and follow the code below.

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
rm(list = ls())
# setwd and read in your datasets:
setwd("~/Desktop/Taxonomic Sequencing/amplicon_bioinformatics/tax_pipe_Mar20/")
dd <- "~/Desktop/Taxonomic Sequencing/amplicon_bioinformatics/package_deal/"

# trait and tax mapping:
source(paste0(dd,"traitmapper_Ramond.R"))
source(paste0(dd,"taxmapper.R"))
source(paste0(dd,"pr2_tax_miner.R"))
source(paste0(dd,"analyze_traitmap_byTrait.R"))

# helper fcns:
source(paste0(dd,"LCA2df.R"))
source(paste0(dd,"idtax2df_pr2.R"))
source(paste0(dd,"bayestax2df.R"))
source(paste0(dd,"idtax2df_silva.R"))
source(paste0(dd,"find_asvs_by_name.R"))

# ensemble tax algrithms:
source(paste0(dd,"consensus_tax_LCAlike.R"))
source(paste0(dd,"consensus_tax_bestRez.R"))
source(paste0(dd,"compare_taxrez.R"))
bayes.pr2 <- readRDS("initial_tax_tabs/bayes_pr2_0boot_Mar20.rds")
idtax.pr2 <- readRDS("initial_tax_tabs/idtax_pr2_0boot_Mar20.rds")
bayes.silva <- readRDS("initial_tax_tabs/bayes_silva_0boot_Mar20.rds")
idtax.silva <- readRDS("initial_tax_tabs/idtax_silva_0boot_Mar20.rds")
lca.pr2 <- readRDS("initial_tax_tabs/LCA_pr2_rawdf_Mar20.rds")
lca.silva <- readRDS("initial_tax_tabs/LCA_silva_rawdf_Mar20.rds")
```

### Arranging and formating our taxonomy tables for running the algorithm:

The data we’re using was pulled slightly haphazardly, so here we’ll use
some bootstrapping estimates to NA-out low-confidence assignments,
reformat our taxonomy tables as dataframes, and sort them alphabetically
by ASV sequences so that the order of rows/ASVs is the same across all
taxonomy
tables.

``` r
# here's the rubric for aligning ASV numbers and sequences across datasets:
library("DECIPHER")
```

    ## Loading required package: Biostrings

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: RSQLite

``` r
setwd("~/Desktop/Taxonomic Sequencing/amplicon_bioinformatics/tax_pipe_Mar20/")
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
```

You can run this for a sanity check:

``` r
# check that they're all in the same order
identical(bayes.pr2$ASV, bayes.silva$ASV)
identical(bayes.pr2$ASV, idtax.pr2$ASV)
identical(idtax.pr2$ASV, idtax.silva$ASV)  
identical(idtax.silva$ASV, lca.pr2$ASV)
identical(lca.pr2$ASV, lca.silva$ASV)
```

…and this to see what the data sets look like. These data sets are
available in the test-data directory.

``` r
head(bayes.pr2)
head(bayes.silva)
head(idtax.pr2)
head(idtax.silva)
head(lca.pr2)
head(lca.silva)
```

The algorithm expects that the ASV’s in each taxonomy data frame are all
in the same order. Before inputting our tables in the algorithm, we will
standardize our `bayes.silva`, `idtax.silva`, and `lca.silva` taxonomies
to pr2 using the taxonomy mapping function as well as sort them by svN
and
ASV’s.

``` r
source("~/Desktop/Taxonomic Sequencing/amplicon_bioinformatics/taxonomy_pipeline/tax_table_mapping/taxmapper.R")

synonym.filepath <- "~/Desktop/Taxonomic Sequencing/amplicon_bioinformatics/taxonomy_pipeline/tax_table_mapping/tax_synonyms_FINAL.csv"

# Bacteria and Archaea doesn't exist in pr2
nonexistent <- c('Bacteria', 'Archaea')

pr2 <- read.csv("~/Desktop/Taxonomic Sequencing/amplicon_bioinformatics/taxonomy_pipeline/tax_table_mapping/pr2_all_tax.csv")
pr2 <- pr2[,-1]

bayes.silva.2.pr2 <- taxmapper(taxin=bayes.silva, tax2map2=pr2, 
                               exceptions=nonexistent,
                               synonym.file=synonym.filepath)

idtax.silva.2.pr2 <- taxmapper(taxin=idtax.silva, tax2map2=pr2, 
                               exceptions=nonexistent,
                               synonym.file=synonym.filepath)
lca.silva.2.pr2 <- taxmapper(taxin=lca.silva, tax2map2=pr2, 
                               exceptions=nonexistent,
                               synonym.file=synonym.filepath)

bayes <- bayes.silva.2.pr2[[3]]
idtax <- idtax.silva.2.pr2[[3]]
lca <- lca.silva.2.pr2[[3]]

s.bayes <- bayes[sort(as.character(bayes$svN), index.return=TRUE)$ix, ]
s.idtax <- idtax[sort(as.character(idtax$svN), index.return=TRUE)$ix, ]
s.lca <- lca[sort(as.character(lca$svN), index.return=TRUE)$ix, ]
```

### First run

Our data should be good to go, so let’s run the algorithm. We have to
specify names for our taxonmy tables in a character vector. We won’t
specify rank names b/c the defaults work for us. There is an option for
the algorithm to count the term NA as part of the majority computation.
In this case, we’ll count it.

First we’ll tell R where to find the algorithm and load it into our
session. In our first run, we’ll specify no tie-breakers, so that the
algorithm will only assign consensus taxonomies to ASVs where a single
taxonomy table has the most common taxonomy of our 3
datasets.

``` r
source("~/Desktop/Taxonomic Sequencing/amplicon_bioinformatics/taxonomy_pipeline/consensus_taxonomies/consensus_tax_mostCom3.R")

tblnam <- c("bayes-pr2", "idtax-pr2", "lca-pr2")
test1 <- consensus_tax_mostCom3(s.bayes, s.idtax, s.lca,
                        tablenames = tblnam, count.na=TRUE,
                        tiebreakz = "none")
head(test1)
```

    ##       svN
    ## 1     sv1
    ## 2    sv10
    ## 3   sv100
    ## 4  sv1000
    ## 5 sv10000
    ## 6 sv10001
    ##                                                                                                                                   ASV
    ## 1     GCTACTACCGATTGAACATTTTAGTGAGGTCCTCGGACTGTGAGCCAGGCGGGTCGCCCTGCCTGGTCTACGGGAAGACGACCAAACTGTAGTGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC
    ## 2     GCACCTACCGATTGAATGGTCCGGTGAAGCCTCGGGATTGTGATTAGTTTCCTTTATTGGAAGGTAGTTATGAGAACCTGTCTAAACCTTATCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCT
    ## 3     GCTCCTACCGATTGAGTGGTCCGGTGAATAATTCGGACTGGTGCCGATTTCGGTTCTCCGAGTTCGGCGCTGGGAAGTCTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 4 GCTCCTACCGATTGAATGGTCCGGTGAAGTGTTCGGATCGTGGCGACGTGGGCGGTTCGCTGCCTGCGACGTCGCGAGAAGTCCACTGAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 5                                  GTTTCTTCCGACTAATGCTTTATGTGAGTGTCACGGATTTTAAAGAAGTGCTGTGAACATTGAGTATCGGAGGAAGAAAAAGTCGTAACAAGGTTATC
    ## 6      GCTGCTACCGATTGAGTGTCCTGGTGAATTATTTGGACCGGCAGTAATTCGAGTTTCTCGATTTACAGCTGGAAAATCTTGTAAACCCTGACACTTAGAGGAAGCAGAAGTCGTAACAAGGTTTCC
    ##     kingdom    supergroup       division           class             order
    ## 1 Eukaryota  Opisthokonta        Metazoa      Arthropoda         Crustacea
    ## 2 Eukaryota Stramenopiles     Ochrophyta Bacillariophyta Bacillariophyta_X
    ## 3 Eukaryota     Alveolata Dinoflagellata     Syndiniales              <NA>
    ## 4 Eukaryota          <NA>           <NA>            <NA>              <NA>
    ## 5      <NA>          <NA>           <NA>            <NA>              <NA>
    ## 6 Eukaryota     Alveolata Dinoflagellata            <NA>              <NA>
    ##           family            genus species
    ## 1    Maxillopoda             <NA>    <NA>
    ## 2 Raphid-pennate Pseudo-nitzschia    <NA>
    ## 3           <NA>             <NA>    <NA>
    ## 4           <NA>             <NA>    <NA>
    ## 5           <NA>             <NA>    <NA>
    ## 6           <NA>             <NA>    <NA>

The output is just a consensus taxonomy table.

``` r
library(readxl)
consensus.res <- read_excel("~/Desktop/Taxonomic Sequencing/amplicon_bioinformatics/taxonomy_pipeline/test_data/bayes_idtax_lca_pr2_mapped.xlsx")
```

    ## New names:
    ## * `` -> ...3
    ## * `` -> ...12

``` r
make.true.NA <- function(x) if(is.character(x)||is.factor(x)){
    is.na(x) <- x=="NA"; x} else {
      x}

consensus.res[] <- lapply(consensus.res, make.true.NA)

colnames(consensus.res)[12] <- "table"

head(consensus.res)
```

    ## # A tibble: 6 x 12
    ##   svN   ASV   ...3  kingdom supergroup division class order family genus species
    ##   <chr> <chr> <lgl> <chr>   <chr>      <chr>    <chr> <chr> <chr>  <chr> <chr>  
    ## 1 sv24… GCTT… NA    Archaea <NA>       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>   
    ## 2 <NA>  <NA>  NA    Archaea <NA>       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>   
    ## 3 <NA>  <NA>  NA    Archaea <NA>       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>   
    ## 4 <NA>  <NA>  NA    Archaea <NA>       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>   
    ## 5 <NA>  <NA>  NA    <NA>    <NA>       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>   
    ## 6 sv14… ACAC… NA    Bacter… <NA>       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>   
    ## # … with 1 more variable: table <chr>

In order to show different scenarios of our algorithm resolving the
majorities, we’ve handpicked certain rows to build a subset data frame
to showcase our
algorithm.

``` r
svs <- c("sv24097", "sv14101", "sv15679", "sv19788", "sv22921", "sv9968", "sv19154", "sv17559", "sv19589", "sv17897", "sv17874", "sv104", "sv5216", "sv5759", "sv6408", "sv5605", "sv12184", "sv16779", "sv15370", "sv20943", "sv23365", "sv1509", "sv4792", "sv22976", "sv20635", "sv17488", "sv17511", "sv22291", "sv8232", "sv15734", "sv5756", "sv14822", "sv21647", "sv4356", "sv9792", "sv5337", "sv11054", "sv5272")

mini.b <- s.bayes[which(s.bayes$svN %in% svs), ]
mini.i <- s.idtax[which(s.idtax$svN %in% svs), ]
mini.l <- s.lca[which(s.lca$svN %in% svs), ]

c.test <- consensus_tax_mostCom3(mini.b, mini.i, mini.l, 
                                 tablenames = c("bayes", "idtax", "lca"), 
                                 ranknamez = c("kingdom", "supergroup", "division", "class", "order", "family", "genus", "species"), 
                                 tiebreakz = "none", count.na = TRUE)
head(c.test)
```

    ##       svN
    ## 1   sv104
    ## 2 sv11054
    ## 3 sv12184
    ## 4 sv14101
    ## 5 sv14822
    ## 6  sv1509
    ##                                                                                                                                    ASV
    ## 1                               GTTTCTTCCGACTGATACTCATTGTGAGTTGCAAGGACCTGTAATGGGAAATTGCTGCAAATGCTTTGTATTGGAGGAAGAAAAAGTCGTAACAAGGTTATC
    ## 2    GCTACTACCGATTGAATGGCTCAGTGAGGCGTTCGGACTGGCCCAGGGAGGTCGGCAACGACCACCCAGGGCCGGAAAGTTCGTCAAACTTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGTCTCC
    ## 3         GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 4                             ACACCACGAGAGTTGGCCGTGCCCGATGTCGTGACTCCAACCCTTCGGGAGGGGAGCGCCTACGGCAAGGTCGGCGATTGGGGTGAAGTCGTAACAAGGTAGCC
    ## 5        GTTGCTACCGATTGATCTGCTGGTAGAGATTGGCCGACTCGGGGATTTAGTGGCAACGCTACTTTTCTGGGGAAACCAATCAATATCGTCAGGTTAGAGGAAGCAAAAGTCGTAACAAGGTTGCT
    ## 6 GCTACTACTGATTGAATTATTTAGTGAGGTCTCCGGACGTGATCACTGTGACGCTTCTAGTGTTACGGTTGTTTCGCAAAAGTTGACCGAACTTGATTATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC
    ##     kingdom   supergroup       division       class          order
    ## 1 Eukaryota         <NA>           <NA>        <NA>           <NA>
    ## 2 Eukaryota Opisthokonta          Fungi  Ascomycota Pezizomycotina
    ## 3 Eukaryota    Alveolata Dinoflagellata Dinophyceae           <NA>
    ## 4  Bacteria         <NA>           <NA>        <NA>           <NA>
    ## 5 Eukaryota         <NA>           <NA>        <NA>           <NA>
    ## 6 Eukaryota Opisthokonta        Metazoa  Arthropoda       Hexapoda
    ##            family   genus species
    ## 1            <NA>    <NA>    <NA>
    ## 2 Dothideomycetes    <NA>    <NA>
    ## 3            <NA>    <NA>    <NA>
    ## 4            <NA>    <NA>    <NA>
    ## 5            <NA>    <NA>    <NA>
    ## 6         Insecta Diptera    <NA>

We’ve saved the result in a csv file called
`bayes_idtax_lca_pr2_mapped.csv` that we will store into a data frame
called `consensus.res`. We chunked out each ASV’s taxonomic assignments
from each of our 3 input taxonomy tables, as well as our consensus
taxonomy table, and put them in consecutive rows. This way, we can look
at chunks of assignments for the same ASV to `consensus.res` the
algorithm and ensure it’s doing what we think.

Here, we’ll pop out a few good examples.

``` r
cbind(rbind(as.data.frame(consensus.res[1, c(1:2, 4:11)]),
      cbind(consensus.res[1, 1:2], consensus.res[2, 4:11]),
      cbind(consensus.res[1, 1:2], consensus.res[3, 4:11]),
      consensus_tax_mostCom3(as.data.frame(consensus.res[1, c(1:2, 4:11)]), 
                       cbind(consensus.res[1, 1:2], consensus.res[2, 4:11]), 
                       cbind(consensus.res[1, 1:2], consensus.res[3, 4:11]),
                       tablenames=c('bayes', 'idtax', 'lca'), 
                       tiebreakz='none', count.na=TRUE)),consensus.res[1:4, 12])
```

    ##       svN
    ## 1 sv24097
    ## 2 sv24097
    ## 3 sv24097
    ## 4 sv24097
    ##                                                                                             ASV
    ## 1 GCTTCACCCGAGTTGGGTTTGGGGGAGATAGTGTCTTATTGACATTATCGAACCTAGGTTCGACAAGGGGGGAGAAGTCGTAACAAGGTGGCC
    ## 2 GCTTCACCCGAGTTGGGTTTGGGGGAGATAGTGTCTTATTGACATTATCGAACCTAGGTTCGACAAGGGGGGAGAAGTCGTAACAAGGTGGCC
    ## 3 GCTTCACCCGAGTTGGGTTTGGGGGAGATAGTGTCTTATTGACATTATCGAACCTAGGTTCGACAAGGGGGGAGAAGTCGTAACAAGGTGGCC
    ## 4 GCTTCACCCGAGTTGGGTTTGGGGGAGATAGTGTCTTATTGACATTATCGAACCTAGGTTCGACAAGGGGGGAGAAGTCGTAACAAGGTGGCC
    ##   kingdom supergroup division class order family genus species  table
    ## 1 Archaea       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>  Bayes
    ## 2 Archaea       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>  Idtax
    ## 3 Archaea       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>    LCA
    ## 4 Archaea       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA> RESULT

The output shows the taxonomic assignments for an arbitrary ASV in 3 of
our initial taxonomy tables, as well as our consensus table. The last
column shows the taxonomy table of each row for this arbitrary ASV.
RESULT deontes our consensus taxonomy. We can see that all of the
taxonomies are the same, so it’s obvious that the resulting consensus
taxonomy should be the same as well.

Now let’s look at one where some contain different taxonomies.

``` r
cbind(rbind(as.data.frame(consensus.res[11, c(1:2, 4:11)]),
      cbind(consensus.res[11, 1:2], consensus.res[12, 4:11]),
      cbind(consensus.res[11, 1:2], consensus.res[13, 4:11]),
      consensus_tax_mostCom3(as.data.frame(consensus.res[11, c(1:2, 4:11)]), 
                       cbind(consensus.res[11, 1:2], consensus.res[12, 4:11]), 
                       cbind(consensus.res[11, 1:2], consensus.res[13, 4:11]),
                       tablenames=c('bayes', 'idtax', 'lca'), 
                       tiebreakz='none', count.na=TRUE)),consensus.res[11:14, 12])
```

    ##       svN
    ## 1 sv15679
    ## 2 sv15679
    ## 3 sv15679
    ## 4 sv15679
    ##                                                                                                                                                    ASV
    ## 1 GCTACTACCGATTGAATGGTCCGGTGAAATCCTCGGAGCCGTGGCCTCTACGCAATCCGGGCAACCGGGTTGTGAGGTCTCCCCTTTTGGCGGCGAAGTCGATTGAACCTTACCATTTAGAGGAAGGAGAAGTCGTAACAAGGTCTCC
    ## 2 GCTACTACCGATTGAATGGTCCGGTGAAATCCTCGGAGCCGTGGCCTCTACGCAATCCGGGCAACCGGGTTGTGAGGTCTCCCCTTTTGGCGGCGAAGTCGATTGAACCTTACCATTTAGAGGAAGGAGAAGTCGTAACAAGGTCTCC
    ## 3 GCTACTACCGATTGAATGGTCCGGTGAAATCCTCGGAGCCGTGGCCTCTACGCAATCCGGGCAACCGGGTTGTGAGGTCTCCCCTTTTGGCGGCGAAGTCGATTGAACCTTACCATTTAGAGGAAGGAGAAGTCGTAACAAGGTCTCC
    ## 4 GCTACTACCGATTGAATGGTCCGGTGAAATCCTCGGAGCCGTGGCCTCTACGCAATCCGGGCAACCGGGTTGTGAGGTCTCCCCTTTTGGCGGCGAAGTCGATTGAACCTTACCATTTAGAGGAAGGAGAAGTCGTAACAAGGTCTCC
    ##     kingdom supergroup division    class         order          family
    ## 1 Eukaryota  Amoebozoa   Lobosa Lobosa_X Centramoebida Acanthamoebidae
    ## 2 Eukaryota  Amoebozoa   Lobosa Lobosa_X Centramoebida Acanthamoebidae
    ## 3 Eukaryota       <NA>     <NA>     <NA>          <NA>            <NA>
    ## 4 Eukaryota  Amoebozoa   Lobosa Lobosa_X Centramoebida Acanthamoebidae
    ##          genus species  table
    ## 1 Acanthamoeba    <NA>  Bayes
    ## 2 Acanthamoeba    <NA>  Idtax
    ## 3         <NA>    <NA>    LCA
    ## 4 Acanthamoeba    <NA> RESULT

Here, we can see that the taxonomies in bayes and idtax are the majority
in each column compared to the NA taxonomy in LCA. Therefore, the
resulting consensus taxonomy would be the taxonomies from bayes and
idtax.

``` r
cbind(rbind(as.data.frame(consensus.res[81, c(1:2, 4:11)]),
      cbind(consensus.res[81, 1:2], consensus.res[82, 4:11]),
      cbind(consensus.res[81, 1:2], consensus.res[83, 4:11]),
      consensus_tax_mostCom3(as.data.frame(consensus.res[81, c(1:2, 4:11)]), 
                       cbind(consensus.res[81, 1:2], consensus.res[82, 4:11]), 
                       cbind(consensus.res[81, 1:2], consensus.res[83, 4:11]),
                       tablenames=c('bayes', 'idtax', 'lca'), 
                       tiebreakz='none', count.na=TRUE)),consensus.res[81:84, 12])
```

    ##       svN
    ## 1 sv12184
    ## 2 sv12184
    ## 3 sv12184
    ## 4 sv12184
    ##                                                                                                                            ASV
    ## 1 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 2 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 3 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 4 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ##     kingdom supergroup       division            class               order
    ## 1 Eukaryota  Alveolata    Apicomplexa Gregarinomorphea Cryptogregarinorida
    ## 2 Eukaryota  Alveolata Dinoflagellata      Dinophyceae       Gonyaulacales
    ## 3 Eukaryota  Alveolata Dinoflagellata      Dinophyceae                <NA>
    ## 4 Eukaryota  Alveolata Dinoflagellata      Dinophyceae                <NA>
    ##              family           genus species  table
    ## 1 Cryptosporidiidae Cryptosporidium    <NA>  Bayes
    ## 2   Goniodomataceae     Alexandrium    <NA>  Idtax
    ## 3              <NA>            <NA>    <NA>    LCA
    ## 4              <NA>            <NA>    <NA> RESULT

We can see that in the Order, Family, and Genus rankings, there is no
majority as there is a three way tie between the taxonomy tables. Since
we didn’t specify a tiebreaker preference, the algorithm automatically
assigns NA to that column as it wasn’t able to conclude to a consensus.

### Second run

Here we will introduce how our algorithm takes care of tie breakers set
by the user input. Let’s say we want to prioritize the bayes taxonomy
table as the tie breaker. This is specified by introducing a list of
pairs of (table name, taxonomy) in the order of priority. This means
that the first element in the list will be prioritized the most while
the last element in thel ist will be prioritized the least.

``` r
cbind(rbind(as.data.frame(consensus.res[81, c(1:2, 4:11)]),
      cbind(consensus.res[81, 1:2], consensus.res[82, 4:11]),
      cbind(consensus.res[81, 1:2], consensus.res[83, 4:11]),
      consensus_tax_mostCom3(as.data.frame(consensus.res[81, c(1:2, 4:11)]), 
                       cbind(consensus.res[81, 1:2], consensus.res[82, 4:11]), 
                       cbind(consensus.res[81, 1:2], consensus.res[83, 4:11]),
                       tablenames=c('bayes', 'idtax', 'lca'), 
                       tiebreakz=list(c('bayes', NA)), count.na=TRUE)),consensus.res[81:84, 12])
```

    ##       svN
    ## 1 sv12184
    ## 2 sv12184
    ## 3 sv12184
    ## 4 sv12184
    ##                                                                                                                            ASV
    ## 1 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 2 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 3 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 4 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ##     kingdom supergroup       division            class               order
    ## 1 Eukaryota  Alveolata    Apicomplexa Gregarinomorphea Cryptogregarinorida
    ## 2 Eukaryota  Alveolata Dinoflagellata      Dinophyceae       Gonyaulacales
    ## 3 Eukaryota  Alveolata Dinoflagellata      Dinophyceae                <NA>
    ## 4 Eukaryota  Alveolata Dinoflagellata      Dinophyceae Cryptogregarinorida
    ##              family           genus species  table
    ## 1 Cryptosporidiidae Cryptosporidium    <NA>  Bayes
    ## 2   Goniodomataceae     Alexandrium    <NA>  Idtax
    ## 3              <NA>            <NA>    <NA>    LCA
    ## 4 Cryptosporidiidae Cryptosporidium    <NA> RESULT

We can see that whenever there is a three way tie, the taxonomy from the
bayes table is assigned as we specified for it to be the tie breaker.

The same thing can be applied if we want to prioritize the idtax table.

``` r
cbind(rbind(as.data.frame(consensus.res[81, c(1:2, 4:11)]),
      cbind(consensus.res[81, 1:2], consensus.res[82, 4:11]),
      cbind(consensus.res[81, 1:2], consensus.res[83, 4:11]),
      consensus_tax_mostCom3(as.data.frame(consensus.res[81, c(1:2, 4:11)]), 
                       cbind(consensus.res[81, 1:2], consensus.res[82, 4:11]), 
                       cbind(consensus.res[81, 1:2], consensus.res[83, 4:11]),
                       tablenames=c('bayes', 'idtax', 'lca'), 
                       tiebreakz=list(c('idtax', NA)), count.na=TRUE)),consensus.res[81:84, 12])
```

    ##       svN
    ## 1 sv12184
    ## 2 sv12184
    ## 3 sv12184
    ## 4 sv12184
    ##                                                                                                                            ASV
    ## 1 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 2 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 3 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 4 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ##     kingdom supergroup       division            class               order
    ## 1 Eukaryota  Alveolata    Apicomplexa Gregarinomorphea Cryptogregarinorida
    ## 2 Eukaryota  Alveolata Dinoflagellata      Dinophyceae       Gonyaulacales
    ## 3 Eukaryota  Alveolata Dinoflagellata      Dinophyceae                <NA>
    ## 4 Eukaryota  Alveolata Dinoflagellata      Dinophyceae       Gonyaulacales
    ##              family           genus species  table
    ## 1 Cryptosporidiidae Cryptosporidium    <NA>  Bayes
    ## 2   Goniodomataceae     Alexandrium    <NA>  Idtax
    ## 3              <NA>            <NA>    <NA>    LCA
    ## 4   Goniodomataceae     Alexandrium    <NA> RESULT

We can prioritize specific pairs of taxonomy and its correspnoding
taxonomy table by replacing the NA with taxonomies we are interested in.
If we want to prioritize the name “NA”, we would put in “na”.

For example, we’ll prioritize the bayes table whenever we encounter the
taxonomy Cryptosporidiidae. One thing to note that the tiebreaking is
done by exact name matching which means that it is also case sensitive.

``` r
cbind(rbind(as.data.frame(consensus.res[81, c(1:2, 4:11)]),
      cbind(consensus.res[81, 1:2], consensus.res[82, 4:11]),
      cbind(consensus.res[81, 1:2], consensus.res[83, 4:11]),
      consensus_tax_mostCom3(as.data.frame(consensus.res[81, c(1:2, 4:11)]), 
                       cbind(consensus.res[81, 1:2], consensus.res[82, 4:11]), 
                       cbind(consensus.res[81, 1:2], consensus.res[83, 4:11]),
                       tablenames=c('bayes', 'idtax', 'lca'), 
                       tiebreakz=list(c('bayes', 'Cryptosporidiidae')), count.na=TRUE)),consensus.res[81:84, 12])
```

    ##       svN
    ## 1 sv12184
    ## 2 sv12184
    ## 3 sv12184
    ## 4 sv12184
    ##                                                                                                                            ASV
    ## 1 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 2 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 3 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 4 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ##     kingdom supergroup       division            class               order
    ## 1 Eukaryota  Alveolata    Apicomplexa Gregarinomorphea Cryptogregarinorida
    ## 2 Eukaryota  Alveolata Dinoflagellata      Dinophyceae       Gonyaulacales
    ## 3 Eukaryota  Alveolata Dinoflagellata      Dinophyceae                <NA>
    ## 4 Eukaryota  Alveolata Dinoflagellata      Dinophyceae                <NA>
    ##              family           genus species  table
    ## 1 Cryptosporidiidae Cryptosporidium    <NA>  Bayes
    ## 2   Goniodomataceae     Alexandrium    <NA>  Idtax
    ## 3              <NA>            <NA>    <NA>    LCA
    ## 4 Cryptosporidiidae            <NA>    <NA> RESULT

We can see here that in the genus, family, and order ranks, we encourted
three way ties. However, in the family rank, it contained the taxonomy
Cryptosporidiidae in the ties, so we prioritized it to be the
tiebreaker. There was no tiebreaking specified for the others, so by
default, it resorted to NA.

### Third run

Here, we’ll introduce more complex tiebreaking inputs that the algorithm
can execute. Remember that the priorities are expressed based on the
order of the pairs of table name and taxonomy in the list. If a certain
taxonomy table is prioritized, its pair would be (table name, NA). If
the specific taxonomy name “NA” wants to be prioritized, then its pair
would be (table name, “na”).

For example, let’s build on what we had earlier and combine our tie
breaking priorities. We will always prioritize the specific pair of
(‘bayes’, ‘Cryptosporidiidae’) and then the idtax table itself.

``` r
cbind(rbind(as.data.frame(consensus.res[81, c(1:2, 4:11)]),
      cbind(consensus.res[81, 1:2], consensus.res[82, 4:11]),
      cbind(consensus.res[81, 1:2], consensus.res[83, 4:11]),
      consensus_tax_mostCom3(as.data.frame(consensus.res[81, c(1:2, 4:11)]), 
                       cbind(consensus.res[81, 1:2], consensus.res[82, 4:11]), 
                       cbind(consensus.res[81, 1:2], consensus.res[83, 4:11]),
                       tablenames=c('bayes', 'idtax', 'lca'), 
                       tiebreakz=list(c('bayes', 'Cryptosporidiidae'), c('idtax', NA)),
                       count.na=TRUE)),consensus.res[81:84, 12])
```

    ##       svN
    ## 1 sv12184
    ## 2 sv12184
    ## 3 sv12184
    ## 4 sv12184
    ##                                                                                                                            ASV
    ## 1 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 2 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 3 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 4 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ##     kingdom supergroup       division            class               order
    ## 1 Eukaryota  Alveolata    Apicomplexa Gregarinomorphea Cryptogregarinorida
    ## 2 Eukaryota  Alveolata Dinoflagellata      Dinophyceae       Gonyaulacales
    ## 3 Eukaryota  Alveolata Dinoflagellata      Dinophyceae                <NA>
    ## 4 Eukaryota  Alveolata Dinoflagellata      Dinophyceae       Gonyaulacales
    ##              family           genus species  table
    ## 1 Cryptosporidiidae Cryptosporidium    <NA>  Bayes
    ## 2   Goniodomataceae     Alexandrium    <NA>  Idtax
    ## 3              <NA>            <NA>    <NA>    LCA
    ## 4 Cryptosporidiidae     Alexandrium    <NA> RESULT

From eariler runs, we noticed that the ranks: order, family, and genus
encountered three way ties. The genus and order ranks got resolved by
prioritizing the idtax table. On the other hand, in the family rank,
since one of the ties was Cryptosporidiidae from the bayes table, it got
prioritized over the idtax table due to the order of pairs we specified
in the list.

Overall, a good tiebreak input would contain specific pairs of table and
taxonomy and always contain a table to prioritize over in the end to
avoid NA assignment.

### Fourth Run

Finally, weights can be introduced into the algorithm to mitigate any
ties. This is also another way to introduce bias toward taxonomy table
when resolving a consensus taxonomy. Here, we’ll introduce a greater
weight on the bayes taxonomy table where it will be inputted twice in
the pool of taxonomies to resolve a consensus to.

``` r
consensus_tax_mostCom3(as.data.frame(consensus.res[81, c(1:2, 4:11)]), 
                       cbind(consensus.res[81, 1:2], consensus.res[82, 4:11]), 
                       cbind(consensus.res[81, 1:2], consensus.res[83, 4:11]),
                       tablenames=c('bayes', 'idtax', 'lca'), 
                       count.na=TRUE, weights=c(2,1,1)) 
```

    ##       svN
    ## 1 sv12184
    ##                                                                                                                            ASV
    ## 1 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ##     kingdom supergroup division class               order            family
    ## 1 Eukaryota  Alveolata     <NA>  <NA> Cryptogregarinorida Cryptosporidiidae
    ##             genus species
    ## 1 Cryptosporidium    <NA>

``` r
cbind(rbind(as.data.frame(consensus.res[81, c(1:2, 4:11)]),
      cbind(consensus.res[81, 1:2], consensus.res[82, 4:11]),
      cbind(consensus.res[81, 1:2], consensus.res[83, 4:11]),
      consensus_tax_mostCom3(as.data.frame(consensus.res[81, c(1:2, 4:11)]), 
                       cbind(consensus.res[81, 1:2], consensus.res[82, 4:11]), 
                       cbind(consensus.res[81, 1:2], consensus.res[83, 4:11]),
                       tablenames=c('bayes', 'idtax', 'lca'), 
                       count.na=TRUE, weights=c(2,1,1))),consensus.res[81:84, 12])
```

    ##       svN
    ## 1 sv12184
    ## 2 sv12184
    ## 3 sv12184
    ## 4 sv12184
    ##                                                                                                                            ASV
    ## 1 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 2 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 3 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 4 GCTCCTACCGATTGAATGATCCGGTGAATAATTCGGACTGGGAAATTTTTAGTTTCTATTCTTTTCACGGGAAGTTTAATAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ##     kingdom supergroup       division            class               order
    ## 1 Eukaryota  Alveolata    Apicomplexa Gregarinomorphea Cryptogregarinorida
    ## 2 Eukaryota  Alveolata Dinoflagellata      Dinophyceae       Gonyaulacales
    ## 3 Eukaryota  Alveolata Dinoflagellata      Dinophyceae                <NA>
    ## 4 Eukaryota  Alveolata           <NA>             <NA> Cryptogregarinorida
    ##              family           genus species  table
    ## 1 Cryptosporidiidae Cryptosporidium    <NA>  Bayes
    ## 2   Goniodomataceae     Alexandrium    <NA>  Idtax
    ## 3              <NA>            <NA>    <NA>    LCA
    ## 4 Cryptosporidiidae Cryptosporidium    <NA> RESULT

Here, we noticed that the resolved taxonomic assignments is just the
bayes table.

### Fifth Run

Lastly, we can make the algorithm not consider NA as a majority
candidate.

``` r
cbind(rbind(as.data.frame(consensus.res[61, c(1:2, 4:11)]),
      cbind(consensus.res[61, 1:2], consensus.res[62, 4:11]),
      cbind(consensus.res[61, 1:2], consensus.res[63, 4:11]),
      consensus_tax_mostCom3(as.data.frame(consensus.res[61, c(1:2, 4:11)]), 
                       cbind(consensus.res[61, 1:2], consensus.res[62, 4:11]), 
                       cbind(consensus.res[61, 1:2], consensus.res[63, 4:11]),
                       tablenames=c('bayes', 'idtax', 'lca'), 
                       count.na=FALSE)),consensus.res[61:64, 12])
```

    ##      svN
    ## 1 sv5216
    ## 2 sv5216
    ## 3 sv5216
    ## 4 sv5216
    ##                                                                                                                            ASV
    ## 1 GCTCCTACCGATTGAGTGATCCGGTGAGTCCTTTAGAGTGTGTTGCAAGTAGTTTCTACCCGCCGCGTGCGAAGTTGAACAAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 2 GCTCCTACCGATTGAGTGATCCGGTGAGTCCTTTAGAGTGTGTTGCAAGTAGTTTCTACCCGCCGCGTGCGAAGTTGAACAAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 3 GCTCCTACCGATTGAGTGATCCGGTGAGTCCTTTAGAGTGTGTTGCAAGTAGTTTCTACCCGCCGCGTGCGAAGTTGAACAAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 4 GCTCCTACCGATTGAGTGATCCGGTGAGTCCTTTAGAGTGTGTTGCAAGTAGTTTCTACCCGCCGCGTGCGAAGTTGAACAAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ##     kingdom    supergroup       division           class             order
    ## 1 Eukaryota     Alveolata    Apicomplexa   Apicomplexa_X          Coccidia
    ## 2 Eukaryota Stramenopiles     Ochrophyta Bacillariophyta Bacillariophyta_X
    ## 3 Eukaryota     Alveolata Dinoflagellata     Syndiniales              <NA>
    ## 4 Eukaryota     Alveolata           <NA>            <NA>              <NA>
    ##                       family         genus species  table
    ## 1                       <NA>          <NA>    <NA>  Bayes
    ## 2 Polar-centric-Mediophyceae Thalassiosira    <NA>  Idtax
    ## 3                       <NA>          <NA>    <NA>    LCA
    ## 4 Polar-centric-Mediophyceae Thalassiosira    <NA> RESULT

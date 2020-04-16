taxmapper\_demo
================
D Catlett & K Son
4/13/2020

## Overview

Here we step through an example use of the taxmapper.R function. Given a
taxonomy table, it maps it to another taxonomy by exact name-matching,
regardless of rank. Requires that all taxonomy tables contain sVN and
ASV’s at the first two columns. Synonyms for certain taxonomies can be
provided to check if exact name matching doesn’t work. To resolve
nonexistent taxonomies from either taxonomy tables, the function assumes
that `exceptions` covers all nonexistent taxonomies in the highest rank,
so it can populate the rest of the ranks with NA regardless.

## Start ’er up:

We’ll clear out our environment, set our wd, and read in taxonomy
tables: The taxonomy tables used here come from implementations of the
RDP Bayesian classifier, the newer idtaxa algorithm, and MEGAN’s LCA
algorithm against both the Silva and pr2 reference databases. Our
amplicon data set is an 18S-V9 tag sequencing project from the coastal
ocean.

You can do this with any taxonomy tables assuming you format them
properly. To follow along with this demo, grab the taxonomy tables in
the “test\_data” directory of this repository and follow the code below.

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

## Arranging and formatting our taxonomy tables:

The data we’re using was pulled slightly haphazardly, so here we’ll use
some bootstrapping estimats to NA-out low-confidence assignments,
reformat our taxonomy tables as data frames, and sort them
alphabetically by ASV sequences so that the order of rows/ASVs is the
same across all taxonomy
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

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

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

## Run taxmapper

Here we’ll do an example run and see what the inputs and outputs look
like. The function iterates through each row of the `taxin` data frame
to map each row to a row from `tax2map2` data frame. The function
returns a list of three elements. The first element is a data frame
including all unique taxonomic assignments and their corresponding
mapped assignments. The second element is a character vector containing
all names within `taxin` that were unable to be mapped. The last element
is a data frame having the ASV seqs supplied in `taxin` with the mapped
taxonomic assignments.

In this example, we’ll map bayes-silva to
bayes-pr2:

``` r
source("~/Desktop/Taxonomic Sequencing/amplicon_bioinformatics/taxonomy_pipeline/tax_table_mapping/taxmapper.R")

synonym.filepath <- "~/Desktop/Taxonomic Sequencing/amplicon_bioinformatics/taxonomy_pipeline/tax_table_mapping/tax_synonyms_FINAL.csv"

# Bacteria and Archaea doesn't exist in pr2
nonexistent <- c('Bacteria', 'Archaea')

pr2 <- read.csv("~/Desktop/Taxonomic Sequencing/amplicon_bioinformatics/taxonomy_pipeline/tax_table_mapping/pr2_all_tax.csv")

bayes.silva.2.pr2 <- taxmapper(taxin=bayes.silva, tax2map2=pr2, 
                               exceptions=nonexistent,
                               synonym.file=synonym.filepath)
```

`bayes.silva.2.pr2` contains a list of the three outputs discussed
earlier. Let’s look through each element, starting with the unique
taxonomic assignments from bayes-silva to bayes-pr2’s corresponding
mapped
    assignments.

``` r
head(bayes.silva.2.pr2[[1]], 10)
```

    ##      Kingdom          Phylum               Class              Order
    ## 1    Archaea   Nanoarchaeota        Nanoarchaeia    Woesearchaeales
    ## 4    Archaea   Nanoarchaeota        Nanoarchaeia    Woesearchaeales
    ## 5    Archaea   Nanoarchaeota        Nanoarchaeia    Woesearchaeales
    ## 6    Archaea   Nanoarchaeota        Nanoarchaeia    Woesearchaeales
    ## 7    Archaea   Nanoarchaeota        Nanoarchaeia    Woesearchaeales
    ## 12  Bacteria  Proteobacteria Gammaproteobacteria     Pasteurellales
    ## 13   Archaea   Euryarchaeota     Methanobacteria Methanobacteriales
    ## 14   Archaea   Nanoarchaeota        Nanoarchaeia    Woesearchaeales
    ## 15   Archaea   Nanoarchaeota        Nanoarchaeia    Woesearchaeales
    ## 75 Eukaryota Chytridiomycota    Chytridiomycetes  Rhizophlyctidales
    ##                 Family         Genus   kingdom   supergroup division
    ## 1      SCGC_AAA286-E23          <NA>   Archaea         <NA>     <NA>
    ## 4               GW2011          AR20   Archaea         <NA>     <NA>
    ## 5                 <NA>          <NA>   Archaea         <NA>     <NA>
    ## 6               GW2011          AR15   Archaea         <NA>     <NA>
    ## 7         CG1-02-57-44          <NA>   Archaea         <NA>     <NA>
    ## 12     Pasteurellaceae Volucribacter  Bacteria         <NA>     <NA>
    ## 13 Methanobacteriaceae          <NA>   Archaea         <NA>     <NA>
    ## 14      SCGC_AAA011-D5          <NA>   Archaea         <NA>     <NA>
    ## 15   GW2011_GWC1_47_15          <NA>   Archaea         <NA>     <NA>
    ## 75  Rhizophlyctidaceae Rhizophlyctis Eukaryota Opisthokonta    Fungi
    ##              class             order           family         genus species
    ## 1             <NA>              <NA>             <NA>          <NA>      NA
    ## 4             <NA>              <NA>             <NA>          <NA>      NA
    ## 5             <NA>              <NA>             <NA>          <NA>      NA
    ## 6             <NA>              <NA>             <NA>          <NA>      NA
    ## 7             <NA>              <NA>             <NA>          <NA>      NA
    ## 12            <NA>              <NA>             <NA>          <NA>      NA
    ## 13            <NA>              <NA>             <NA>          <NA>      NA
    ## 14            <NA>              <NA>             <NA>          <NA>      NA
    ## 15            <NA>              <NA>             <NA>          <NA>      NA
    ## 75 Chytridiomycota Chytridiomycotina Chytridiomycetes Rhizophlyctis      NA

Here we see any instances of Archaea or Bacteria from bayes-silva gets
mapped to Bacteria on the bayes-pr2 side.

Let’s look at the character vector of taxonomies from bayes-silva that
wasn’t able to be
    mapped.

``` r
head(bayes.silva.2.pr2[[2]], 10)
```

    ##  [1] "Rhizophlyctidaceae"    "Rhizophlyctidales"     "Hypotrichomonadea_fa" 
    ##  [4] "Hypotrichomonadea_or"  "Hypotrichomonadea"     "Pinguiochrysidales_fa"
    ##  [7] "Ochrophyta_ph"         "Hildenbrandiophycidae" "Florideophycidae"     
    ## [10] "Tetramitia"

Let’s look at the data frame of ASV’s from bayes-silva mapped to their
corresponding taxonomies from
    bayes-pr2.

``` r
head(bayes.silva.2.pr2[[3]], 10)
```

    ##                                                                                              ASV
    ## 1  AAACCAACCAAGTAAGATCTAGAAAAGGCATTATTTTTTTGGTAATGTCAATTCTTAATTTTGTGAGGGGGGTTAAGTCGTCACAAGGTATCC
    ## 4    AAACCAACCGAGCAGGCTTTAGATGAGCTTCATTGAGTAAATGGATCGAATCTAAGGTCAGTGAGGTGGGTTAAGTCGTCACAAGGTATCC
    ## 5    AAACCAACCGAGCAGGCTTTAGATGAGTTTTATTGAGTAAATAAGTCGAATCTAAGGTCAGTGAGGTGGGTTAAGTCGTCACAAGGTATCC
    ## 6   AAACCAACCGAGCAGGGTCCGGGTAAACCACGCTTTCTTGAGCATGGTGAATCTAGTCTCAGTGAGGTGGGTTAAGTCGTCACAAGGTAACC
    ## 7    AAACCAACCGAGCAGGGTCTGGGTGAGCCTTATTGAGTAAATAAGCCGAACTTGGGCTCAGTGAGGTGGGATAAGTCGTCACAAGGTATCC
    ## 12  AAACCAACCGAGCCTAGCTTGGGTGAGTATTAGTTTTTTAGCTAATCCAAATCTAAGTTAAGTGAGGTGGGTTAAGTCGTCACAAGGTATCC
    ## 13   AAACCAACCGAGCGGGGTTTGGATGAGATATGGCTTTTAGCTATCTCGAATCTAAATTTCGTGAGGCGGGTTAAGTCGTCACAAGGTACCT
    ## 14         AAACCAACCGAGTTATGAAGGGATGAAGCCCTTATTGGGAAAAATCTTTTCATGACAAGGTAGGTTAAGTCGTCACAAGGTATCT
    ## 15         AAACCAACCGAGTTATGAAGGGATGAAGCTCTCATTGAGAAGAATCTTTTCACGACAAGGTAGGTTAAGTCGTCACAAGGTATCT
    ## 75          AAACCAACCGAGTTGTGTTCTGGCGAGGTTTTTCGAAAACGAACCTTAGCACGACAAGGCAGGTTAAGTCGACACAAGGTATCT
    ##      kingdom   supergroup division           class             order
    ## 1    Archaea         <NA>     <NA>            <NA>              <NA>
    ## 4    Archaea         <NA>     <NA>            <NA>              <NA>
    ## 5    Archaea         <NA>     <NA>            <NA>              <NA>
    ## 6    Archaea         <NA>     <NA>            <NA>              <NA>
    ## 7    Archaea         <NA>     <NA>            <NA>              <NA>
    ## 12  Bacteria         <NA>     <NA>            <NA>              <NA>
    ## 13   Archaea         <NA>     <NA>            <NA>              <NA>
    ## 14   Archaea         <NA>     <NA>            <NA>              <NA>
    ## 15   Archaea         <NA>     <NA>            <NA>              <NA>
    ## 75 Eukaryota Opisthokonta    Fungi Chytridiomycota Chytridiomycotina
    ##              family         genus species
    ## 1              <NA>          <NA>      NA
    ## 4              <NA>          <NA>      NA
    ## 5              <NA>          <NA>      NA
    ## 6              <NA>          <NA>      NA
    ## 7              <NA>          <NA>      NA
    ## 12             <NA>          <NA>      NA
    ## 13             <NA>          <NA>      NA
    ## 14             <NA>          <NA>      NA
    ## 15             <NA>          <NA>      NA
    ## 75 Chytridiomycetes Rhizophlyctis      NA

consensus\_tax\_bestRez Demo
================
D Catlett
March 13, 2020

Overview
--------

Here, we step through implementations of the consensus\_tax\_bestRez.R algorithm to demonstrate proper uses and anticipated results. This algorithm is built to utilize information from multiple taxonomy tables (defined here as a table of ASVs and corresponding taxonomic assignments) in order to improve the resolution and/or accuracy of taxonomic annotations of amplicon sequences. It incorporates information from multiple taxonomy tables and determines a "consensus taxonomy" for each ASV in your data set. The algorithm requires that all taxonomy tables follow the same taxonomic naming and ranking conventions, that the order of columns in each taxonomy table follows the taxonomic ranking heirarchy (e.g., Kingdom = taxtable\[,1\], Species = taxtable\[,ncol(taxtable)\]), and that the order of rows (ASVs) in each of the input taxonomy tables is the same. If these rules are not followed, spurious results are likely.

### consensus\_tax\_bestRez Algorithm Description

consensus\_tax\_bestRez generates a consensus taxonomy by using the most highly-resolved taxonomic annotation across all user-specified input taxonomy tables. In other words, it's designed to get the most resolution possible with no regard for whether multiple taxonomies agree on that resolution. Unresolved taxonomic assignments in each taxonomy table should be indicated by NA as the algorithm makes initial assignments by finding the corresponding rows of the input taxonomy tables which have the least number of NA's. There will be occasions when multiple input taxonomy tables are 'tied' for the best taxonomic resolution for a given ASV; thus, **after** assigning a consensus taxonomy to the ASVs where a single taxonomy provides the best resolution, the algorithm uses a series of user-specified rules to break the remaining ties. Rules you may specify are demonstrated below.

### Start 'er up:

We'll clear out our environment, load in the reshape2 package for later, set our wd, and read in taxonomy tables: The taxonomy tables used here come from implementations of the RDP Bayesian classifier, the new idtaxa algorithm, and MEGAN's LCA algorithm against both the Silva and pr2 reference databases. Our amplicon data set is an 18S-V9 tag sequencing project from the coastal ocean.

You can do this with any taxonomy tables assuming you format them properly. To follow along with this demo, grab the taxonomy tables in the "test\_data" directory of this repository and follow the code below.

``` r
rm(list = ls())

# load up reshape2:
.cran_packages <- c("reshape2")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
sapply(c(.cran_packages), require, character.only = TRUE)

# setwd and read in your datasets:
setwd("~/Documents/R/desktop_ampData_processing/connie_taxonomy_stuff_Mar2020/18sV9_amplicon_sequencing/taxonomy_pipeline/demos_and_validation")

idtax.pr2 <- readRDS("~/Documents/R/desktop_ampData_processing/connie_taxonomy_stuff_Mar2020/18sV9_amplicon_sequencing/taxonomy_pipeline/test_data/idtax_0boot_pr2_all18SAug19.rds")
bayes.pr2 <- readRDS("~/Documents/R/desktop_ampData_processing/connie_taxonomy_stuff_Mar2020/18sV9_amplicon_sequencing/taxonomy_pipeline/test_data/bayes_0boot_pr2_all18SAug19.rds")
bayes.silva <- read.csv("~/Documents/R/desktop_ampData_processing/connie_taxonomy_stuff_Mar2020/18sV9_amplicon_sequencing/taxonomy_pipeline/test_data/bayes_silva_60boot_mapped2pr2_all18SAug19.csv",
                        stringsAsFactors = FALSE)
idtax.silva <- read.csv("~/Documents/R/desktop_ampData_processing/connie_taxonomy_stuff_Mar2020/18sV9_amplicon_sequencing/taxonomy_pipeline/test_data/idtax_silva_0boot_mapped2pr2_all18SAug19.csv",
                        stringsAsFactors = FALSE)
lca.pr2 <- read.csv("~/Documents/R/desktop_ampData_processing/connie_taxonomy_stuff_Mar2020/18sV9_amplicon_sequencing/taxonomy_pipeline/test_data/LCA_pr2_mapped2pr2_all18SAug19.csv",
                        stringsAsFactors = FALSE)
lca.silva <- read.csv("~/Documents/R/desktop_ampData_processing/connie_taxonomy_stuff_Mar2020/18sV9_amplicon_sequencing/taxonomy_pipeline/test_data/LCA_silva_mapped2pr2_all18SAug19_Fixed.csv",
                    stringsAsFactors = FALSE)
```

### Arranging and formatting our taxonomy tables for running the algorithm:

The data we're using was pulled slightly haphazardly, so here we'll use some bootstrapping estimates to NA-out low-confidence assignments, reformat our taxonomy tables as dataframes, and sort them alphabetically by ASV sequences so that the order of rows/ASVs is the same across all taxonomy tables.

``` r
# convert tax tables to dataframes as needed and sort by seq's to get the same order..:
conf <- as.data.frame(bayes.pr2$boot, stringsAsFactors = FALSE)
bayes.pr2 <- as.data.frame(bayes.pr2$tax, stringsAsFactors = FALSE)
bayes.pr2[conf < 60] <- NA

source("~/Documents/R/desktop_ampData_processing/connie_taxonomy_stuff_Mar2020/18sV9_amplicon_sequencing/taxonomy_pipeline/helper_fcns/idtax2df.R")
idtax.pr2 <- idtax2df(idtax.pr2, boot = 60)

# sorting each dataframe by DNA sequences:
ii <- sort(rownames(bayes.pr2), index.return = TRUE)
bayes.pr2 <- bayes.pr2[ii$ix,]
idtax.pr2 <- idtax.pr2[ii$ix,]
jj <- sort(bayes.silva$DNASeq, index.return = TRUE)
bayes.silva <- bayes.silva[jj$ix,2:9]
kk <- sort(idtax.silva$Sequence, index.return = TRUE)
idtax.silva <- idtax.silva[kk$ix,3:10]
ll <- sort(lca.silva$Sequence, index.return = TRUE)
lca.silva <- lca.silva[ll$ix,3:10]
mm <- sort(lca.pr2$Sequence, index.return = TRUE)
lca.pr2 <- lca.pr2[mm$ix,3:10]
```

You can run this for a sanity check:

``` r
# compare the sorted sequence arrays to ensure they're all =:
identical(ii$x, jj$x)
identical(jj$x, kk$x)
identical(kk$x, ll$x)
identical(ll$x,mm$x)
```

...and this to see what the data sets look like. These data sets are available in the test-data directory.

``` r
# one more check:
head(bayes.pr2, n = 10)
head(bayes.silva, n = 10)
head(idtax.pr2, n = 10)
head(idtax.silva, n = 10)
head(lca.pr2, n = 10)
head(lca.silva, n = 10)
```

### First run

Our data should be good to go, so let's run the algorithm. We have to specify names for our taxonomy tables in a character vector (*tblnam* below - note they match the order in which we've supplied them to the algorithm). We won't specify rank names b/c the defaults work for us (see the function documentation).

First we'll tell R where to find the algorithm and load it into our session. In our first run, we'll specify no tie-breakers, so that the algorithm will only assign consensus taxonomies to ASVs where a single taxonomy table has the best taxonomic resolution of our 6 datasets.

``` r
source("~/Documents/R/desktop_ampData_processing/connie_taxonomy_stuff_Mar2020/18sV9_amplicon_sequencing/taxonomy_pipeline/consensus_taxonomies/consensus_tax_bestRez.R")

tblnam <- c("bayes-pr2", "idtax-pr2", "lca-pr2", "bayes-silva", "idtax-silva", "lca-silva")
test1 <- consensus_tax_bestRez(bayes.pr2, idtax.pr2, lca.pr2, bayes.silva, idtax.silva, lca.silva,
                        tablenames = tblnam, 
                        tiebreakz = "none")
```

The output is a list with 3 elements. The 1st element is the consensus taxonomy table, the 2nd element is a list of all the input taxonomies, and the third element is a numeric vector containing the row indices of all ASVs that still require tie-breaking. The below code chunk will manipulate the outputs to create a more manageable subset of outputs we can use for qualitative QC to ensure the algorithm works properly.

``` r
itb <- test1[[3]] # indices of rows needing tiebreakerz (none used above)
allt <- test1[[2]] # all input taxonomy tables
# assign same rank names to all input tables:
ranknamez <- c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species") 
ah <- function(x) {
  colnames(x) <- ranknamez
  return(x)
}
subber <- 1:2000 
allt <- lapply(allt, ah) 
alltsub <- lapply(allt, function(x) x[subber,]) # subset a smaller number of rows for QC
eh <- melt(alltsub) # pops it all into a dataframe...
# adds table names to indicate which taxtable each row came from
for (i in 1:length(tblnam)) {
  eh$L1[eh$L1 == i] <- tblnam[i]
}

# add in merged tax to your QC df and write out a csv that you can QC manually...
mt <- test1[[1]]
mt$L1 <- "merged"
QCer <- rbind(eh,mt[subber,],stringsAsFactors = FALSE)
# below rearranges the combined df by ASV rather than by taxtable...
tot <- max(subber)
strt <- min(subber)
yy <- c()
for (i in strt:tot) {
  xx <- seq(from = i, to = tot * (length(tblnam)+1), by = tot)
  yy <- append(yy,xx)
}
QCer <- QCer[yy,]
```

Above we've created a dataframe, *QCer*, that has chunked out each ASV's taxonomic assignments from each of our 6 input taxonomy tables, as well as our consensus taxonomy table, and put them in consecutive rows. This way, we can look at chunks of assignments for the same ASV to QC the algorithm and ensure it's doing what we think. I did this by saving QCer to a csv file and flipping through the outputs, as below:

``` r
write.csv(QCer, file = "~/Documents/R/desktop_ampData_processing/connie_taxonomy_stuff_Mar2020/QCme_6tables_none.csv")
```

Here I'll pop out a few good examples I found from my .csv file to demo the expected output of the *tiebreakerz = "none"*

``` r
ii <- which(rownames(QCer) %in% "1")
QCer[seq(ii,ii+6,1),]
```

    ##         Kingdom Supergroup   Division       Class          Order
    ## 1     Eukaryota       <NA>       <NA>        <NA>           <NA>
    ## 2001       <NA>       <NA>       <NA>        <NA>           <NA>
    ## 4001       <NA>       <NA>       <NA>        <NA>           <NA>
    ## 6001  Eukaryota   Excavata Metamonada Parabasalia Trichomonadida
    ## 8001       <NA>       <NA>       <NA>        <NA>           <NA>
    ## 10001      <NA>       <NA>       <NA>        <NA>           <NA>
    ## 12001 Eukaryota   Excavata Metamonada Parabasalia Trichomonadida
    ##                            Family       Genus Species          L1
    ## 1                            <NA>        <NA>    <NA>   bayes-pr2
    ## 2001                         <NA>        <NA>    <NA>   idtax-pr2
    ## 4001                         <NA>        <NA>    <NA>     lca-pr2
    ## 6001  Trichomitus-Hypotrichomonas Trichomitus    <NA> bayes-silva
    ## 8001                         <NA>        <NA>    <NA> idtax-silva
    ## 10001                        <NA>        <NA>    <NA>   lca-silva
    ## 12001 Trichomitus-Hypotrichomonas Trichomitus    <NA>      merged

This output shows the taxonomic assignments for an arbitrary ASV in 6 of our initial taxonomy tables, as well as our consensus table. The column *L1* shows the taxonomy table of each row for this arbitrary ASV. *merged* denotes our consensus taxonomy. We can see that the *bayes-silva* taxonomy table had the most highly resolved assignment for this ASV, and that assignment was retained in our merged taxonomy. Nice!

Now let's look at one where multiple taxonomies had equivalent resolutions in their assignments:

``` r
ii <- which(rownames(QCer) %in% "143")
QCer[seq(ii,ii+6,1),]
```

    ##         Kingdom Supergroup Division Class Order Family Genus Species
    ## 143   Eukaryota       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ## 2143       <NA>       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ## 4143       <NA>       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ## 6143   Bacteria       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ## 8143   Bacteria       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ## 10143  Bacteria       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ## 12143      <NA>       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ##                L1
    ## 143     bayes-pr2
    ## 2143    idtax-pr2
    ## 4143      lca-pr2
    ## 6143  bayes-silva
    ## 8143  idtax-silva
    ## 10143   lca-silva
    ## 12143      merged

Here we see that the highest taxonomic resolution, which is only to the Kingdom rank for this ASV, was achieved by multiple input taxonomy tables, and we see that our output *merged* taxonomy has not even been assigned at the Kingdom rank. This is actually what we wanted based on our specification of *tiebreakerz = "none"* when we ran the algorithm. Woohoo!!

### Second run

Now we can start trying to break ties. The most simple option is to use *tiebreakerz = "LCAlike"*. This option is inspired by MEGAN's Least Common Ancestor algorithm (see Huson et al, 2007). Here, this specification tells the algorithm to extract the taxonomic assignments that have resulted in a tie for this ASV (so only the most highly resolved assignments), and find the lowest rank at which they are in agreement. The consensus taxonomic assignment will only retain assignments for ranks that are in agreement amongst the taxonomy tables that were most highly resolved for this ASV; all other ranks will be left unassigned.

To see what this looks like in practice, do this:

``` r
# clearing out your environment:
rm(list = setdiff(ls(), c("bayes.pr2","idtax.pr2", "lca.pr2", "bayes.silva", "idtax.silva", "lca.silva", "consensus_tax_bestRez", "tblnam", "subber","itb")))

tblnam <- c("bayes-pr2", "bayes-silva", "idtax-pr2", "lca-pr2", "idtax-silva", "lca-silva")
test2 <- consensus_tax_bestRez(bayes.pr2[itb,], bayes.silva[itb,], idtax.pr2[itb,], lca.pr2[itb,], idtax.silva[itb,], lca.silva[itb,],
                         tablenames = tblnam, 
                         tiebreakz = "LCAlike")
```

Note that we retained the index vector of rows requiring tiebreakers from our first run and only input those rows from our 6 taxonomy tables - that is just to narrow our search for examples of where the LCA-like specification will actually be used.

Re-arrange outputs as above (not shown for brevity), and have a look at the outputs:

``` r
ii <- which(rownames(QCer) %in% "1")
QCer[seq(ii,ii+6,1),]
```

    ##         Kingdom Supergroup Division Class Order Family Genus Species
    ## 1     Eukaryota       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ## 2001   Bacteria       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ## 4001       <NA>       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ## 6001       <NA>       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ## 8001       <NA>       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ## 10001  Bacteria       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ## 12001      <NA>       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ##                L1
    ## 1       bayes-pr2
    ## 2001  bayes-silva
    ## 4001    idtax-pr2
    ## 6001      lca-pr2
    ## 8001  idtax-silva
    ## 10001   lca-silva
    ## 12001      merged

Here we see that our *merged* taxonomy for this particular ASV was left entirely unassigned. This is because amongst the 3 taxonomy tables that were tied with the highest resolution (*bayes-pr2, bayes-silva*, and *lca-silva*), they had no assignments in common (though 2 of them did, the disagreement with the 3rd prevents a true LCA).

Another example:

``` r
ii <- which(rownames(QCer) %in% "143")
QCer[seq(ii,ii+6,1),]
```

    ##         Kingdom    Supergroup        Division Class   Order   Family      Genus
    ## 143   Eukaryota Stramenopiles Stramenopiles_X  MAST MAST-12 MAST-12D MAST-12D_X
    ## 2143  Eukaryota Stramenopiles Stramenopiles_X  MAST  MAST-3  MAST-3C       <NA>
    ## 4143  Eukaryota Stramenopiles Stramenopiles_X  MAST MAST-12 MAST-12D MAST-12D_X
    ## 6143  Eukaryota Stramenopiles            <NA>  <NA>    <NA>     <NA>       <NA>
    ## 8143  Eukaryota Stramenopiles Stramenopiles_X  MAST MAST-12 MAST-12D       <NA>
    ## 10143 Eukaryota          <NA>            <NA>  <NA>    <NA>     <NA>       <NA>
    ## 12143 Eukaryota Stramenopiles Stramenopiles_X  MAST MAST-12 MAST-12D MAST-12D_X
    ##              Species          L1
    ## 143   MAST-12D_X_sp.   bayes-pr2
    ## 2143            <NA> bayes-silva
    ## 4143  MAST-12D_X_sp.   idtax-pr2
    ## 6143            <NA>     lca-pr2
    ## 8143            <NA> idtax-silva
    ## 10143           <NA>   lca-silva
    ## 12143 MAST-12D_X_sp.      merged

Is this what we expected? Yes! The 2 most highly resolved taxonomy tables agree all the way down to their species assignments, so we retain assignments all the way down to the species rank in our *merged* consensus taxonomy. This is not the LCA assignment across all 6 tables, but LCA assignment across those taxonomy tables that were most highly resolved in their assignments for this particular ASV.

One more:

``` r
ii <- which(rownames(QCer) %in% "464")
QCer[seq(ii,ii+6,1),]
```

    ##         Kingdom     Supergroup       Division          Class            Order
    ## 464   Eukaryota           <NA>           <NA>           <NA>             <NA>
    ## 2464  Eukaryota Archaeplastida   Streptophyta Zygnemophyceae Zygnemophyceae_X
    ## 4464       <NA>           <NA>           <NA>           <NA>             <NA>
    ## 6464  Eukaryota      Alveolata Dinoflagellata    Syndiniales             <NA>
    ## 8464       <NA>           <NA>           <NA>           <NA>             <NA>
    ## 10464 Eukaryota      Alveolata Dinoflagellata    Syndiniales    Dino-Group-II
    ## 12464 Eukaryota           <NA>           <NA>           <NA>             <NA>
    ##                      Family       Genus Species          L1
    ## 464                    <NA>        <NA>    <NA>   bayes-pr2
    ## 2464      Zygnemophyceae_XX  Sirogonium    <NA> bayes-silva
    ## 4464                   <NA>        <NA>    <NA>   idtax-pr2
    ## 6464                   <NA>        <NA>    <NA>     lca-pr2
    ## 8464                   <NA>        <NA>    <NA> idtax-silva
    ## 10464 Dino-Group-II-Clade-4 Amoebophrya    <NA>   lca-silva
    ## 12464                  <NA>        <NA>    <NA>      merged

Again, we get what we expected, which is only a Kingdom assignment in our *merged* consensus taxonomy. This is expected because though we have several highly-resolved taxonomy tables for this ASV, the 2 most highly resolved are only in agreement at the kingdom rank. Sickarooni!

### Third run

We can add more complexity to our tie-breaking rules if we know that certain algorithm/database combo's result in more robust taxonomic assignments overall or for particular taxonomic groups (or if we just have an arbitrary preference for any of the above).

We do this by specifying a list for our tiebreakerz argument. The order of the elements in the list dictates the priority of the specified rules (tiebreakz\[\[1\]\] will be utilized preferentially over tiebreakz\[\[2\]\]). 1x2 character vectors should be prescribed to each element of the list. Each element of these character vectors should have one of your tablenames in position 1, and a taxonomic name, or NA in position 2. If a taxonomic name is supplied, the specified table in entry 1 will be given tie-breaking priority only for assignments containing that name. If position 2 is NA, the specified table in entry 1 will be given tie-breaking priority only for all assignments in which it is tied for the best resolution. You can also specify c("LCAlike",NA) to add the LCA-like behavior shown above to your tie-breaking rules. Note that LCAlike will assign all remaining ASVs that require tiebreakers, so placing rules lower than LCA-like will make them useless.

Let's say that we want to give tie-breaking preference to our *idtaxa-silva* followed by *lca-silva* tables, but only when each has assigned *Bacteria*. This is what that would look like (LCAlike is added at the lowest priority for fun):

``` r
rm(list = setdiff(ls(), c("bayes.pr2","idtax.pr2", "lca.pr2", "bayes.silva", "idtax.silva", "lca.silva", "consensus_tax_bestRez", "tblnam", "subber", "itb")))
tb <- list(c("idtax-silva", "Bacteria"), 
           c("lca-silva","Bacteria"), 
           c("LCAlike",NA))
tblnam <- c("bayes-pr2", "bayes-silva", "idtax-pr2", "lca-pr2", "idtax-silva", "lca-silva")
test3 <- consensus_tax_bestRez(bayes.pr2[itb,], bayes.silva[itb,], idtax.pr2[itb,], lca.pr2[itb,], idtax.silva[itb,], lca.silva[itb,],
                        tablenames = tblnam, 
                        tiebreakz = tb)
```

Let's see if it worked:

``` r
ii <- which(rownames(QCer) %in% "1")
QCer[seq(ii,ii+6,1),]
```

    ##         Kingdom Supergroup Division Class Order Family Genus Species
    ## 1     Eukaryota       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ## 5001   Bacteria       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ## 10001      <NA>       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ## 15001      <NA>       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ## 20001      <NA>       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ## 25001  Bacteria       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ## 30001  Bacteria       <NA>     <NA>  <NA>  <NA>   <NA>  <NA>    <NA>
    ##                L1
    ## 1       bayes-pr2
    ## 5001  bayes-silva
    ## 10001   idtax-pr2
    ## 15001     lca-pr2
    ## 20001 idtax-silva
    ## 25001   lca-silva
    ## 30001      merged

Sure did. We see that our *merged* taxonomy has been assigned to Bacteria. In this case, *lca-silva* was one of our most highly-resolved taxonomies while *idtax-silva* was not. Thus, the tie amonst the 3 most highly resolved taxonomies was broken by the *lca-silva Bacteria* assignment.

Another example:

``` r
ii <- which(rownames(QCer) %in% "971")
QCer[seq(ii,ii+6,1),]
```

    ##         Kingdom Supergroup   Division        Class       Order        Family
    ## 971   Eukaryota       <NA>       <NA>         <NA>        <NA>          <NA>
    ## 5971  Eukaryota  Alveolata Ciliophora Spirotrichea Hypotrichia Hypotrichia_X
    ## 10971      <NA>       <NA>       <NA>         <NA>        <NA>          <NA>
    ## 15971      <NA>       <NA>       <NA>         <NA>        <NA>          <NA>
    ## 20971 Eukaryota  Alveolata Ciliophora Spirotrichea Hypotrichia Hypotrichia_X
    ## 25971  Bacteria       <NA>       <NA>         <NA>        <NA>          <NA>
    ## 30971 Eukaryota  Alveolata Ciliophora Spirotrichea Hypotrichia Hypotrichia_X
    ##             Genus Species          L1
    ## 971          <NA>    <NA>   bayes-pr2
    ## 5971  Bergeriella    <NA> bayes-silva
    ## 10971        <NA>    <NA>   idtax-pr2
    ## 15971        <NA>    <NA>     lca-pr2
    ## 20971 Bergeriella    <NA> idtax-silva
    ## 25971        <NA>    <NA>   lca-silva
    ## 30971 Bergeriella    <NA>      merged

This time, we had 2 taxonomy tables tied for the highest resolution, and the tie was broken by the LCAlike specification. Dope.

### Fourth Run

Now we're gonna get more complicated.

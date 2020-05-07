# this script creates "keys" of sample names to subset different groups of samples from 
# your giant 18S re-run

# basically need to look up your old composites from PB, AC, and Tanika and create 
# sample name subsets that you can use to subset for each of those

rm(list = ls())
setwd("~/Documents/R/amplicon_bioinformatics/tax_pipe_Mar20")
old1 <- readRDS("~/Documents/R/desktop_ampData_processing/MBON_all18s_Mar2019/field_analysis/data/ps_object_protistOnly_20190320.rds")
old2 <- readRDS("~/Documents/R/desktop_ampData_processing/MBON_all18s_Mar2019/preProcessing/tax_filtered/time_series_phyloseq_objects/ps_object_protistOnly_20190320.rds")

ac <- readRDS("~/Documents/R/desktop_ampData_processing/ACIDD_18s_Mar2019/field_analysis/data/fields_protistOnly.rds")
tl <- readRDS("~/Documents/R/desktop_ampData_processing/MBON_all18s_Jan2019_analysis/preProcessing/subset4tanika/tax_filtered_tanika_ps_object_20190205.rds")

unique(sample_data(ac)$Sample_type)

pb <- subset_samples(old1, Sample_type %in% c("Field", "PB"))
pbs <- subset_samples(pb, Depth_m == 0)
pbd <- subset_samples(pb, Depth_m %in% c(30, 75, 150, 300))
tl <- subset_samples(tl, Sample_type == "Field")

pbs.sn <- sample_names(pbs)
pbd.sn <- sample_names(pbd)

ac.sn <- sample_names(ac) # controls were already subsetted so these are good
tl.sn <- sample_names(tl)

allem <- c(pbs.sn, pbd.sn, ac.sn, tl.sn)
ml <- length(allem)
tl.sn <- c(tl.sn, rep(NA, times = ml - length(tl.sn)))
ac.sn <- c(ac.sn, rep(NA, times = ml - length(ac.sn)))
pbs.sn <- c(pbs.sn, rep(NA, times = ml - length(pbs.sn)))
pbd.sn <- c(pbd.sn, rep(NA, times = ml - length(pbd.sn)))

subgroups <- data.frame(allfield = allem, pbsurf = pbs.sn, pbdeep = pbd.sn, AC = ac.sn, TL = tl.sn, stringsAsFactors = FALSE)
saveRDS(subgroups, file = "sample_subsets.rds")

### This script combines multiple seuqencing libraries into a single seqtab
### then removes chimeras and saves and output table.
###

# Merge multiple runs (if necessary)
st1 <- readRDS(paste0(pathout_ISV1,"seqtab.rds"))
st2 <- readRDS(paste0(pathout_ISV2,"seqtab.rds"))
st3 <- readRDS(paste0(pathout_ISV3,"seqtab.rds"))
st4 <- readRDS(paste0(pathout_ISV4,"seqtab.rds"))
st5 <- readRDS(paste0(pathout_ISV5,"seqtab.rds"))
st6 <- readRDS(paste0(pathout_ISV6,"seqtab.rds"))
st7 <- readRDS(paste0(pathout_ISV7,"seqtab.rds"))
st8 <- readRDS(paste0(pathout_ISV8,"seqtab.rds"))
st9 <- readRDS(paste0(pathout_ISV9,"seqtab.rds"))
st10 <- readRDS(paste0(pathout_ISV10,"seqtab.rds"))

st.all <- mergeSequenceTables(st1, st2, st3, st4, st5, st6, st7, st8, st9, st10) # merge all runs

seqtab.nochimez <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

Nreads.nochime <- rowSums(seqtab.nochimez)

dir.create(pathout_RC,recursive = TRUE)
saveRDS(seqtab.nochimez, paste0(pathout_RC,"seqtab.rds")) # CHANGE ME to where you want sequence table saved
write.csv(Nreads.nochime,paste0(pathout_RC,"readsChimeraFiltered.csv"))


# UNDER CONSTRUCTION

# this is a shell script that executes the functions I (+Kevin +Connie) written for my ensemble taxonomy pipeline
# also using it as an outline to track where I'm at from start to finish..

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

#### Step 1: Preliminary mapping of individual tax tables onto trait database [could skip]
# make sure boot-strapping thresholds are honored here...
# re-format as dataframes...
# remove non-protists...
# then write a shell of my trait mapping and analysis functions

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


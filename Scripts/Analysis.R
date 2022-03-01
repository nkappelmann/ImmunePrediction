# -------------------------------------------------------
# Analysis-----------------------------------------------
# -------------------------------------------------------


# 1 Preparation------------------------------------------

## Load packages
library("tidyverse")

## Set Slurm directory
setwd("/home/nkappelmann/OPTIMA/ImmunePrediction")


## Load data
load("./Data/OPTIMA_Cytokine_PreprocessedData.RData")


## Index participants used in analysis
# Index participants with minimum of two BDI observations
ids_with_bdi = dat[, paste0("t", 0:7, "_bdi")] %>% is.na() %>% rowSums < 7

# Index participants with genetic data
ids_with_geneticdata = !is.na(dat$BMI_PRS)

# Index participants with transcriptomic data
ids_with_rnadata = !is.na(dat$rna.PC1) # Present in all individuals

# Index participants with all relevant data
ids_for_analysis = ids_with_bdi & ids_with_geneticdata & ids_with_rnadata

## Load cytokine reference
load("./Data/Cytokine_Reference.RData")


## Define üredictor variable sets
prs.vars = colnames(dat)[grepl("PRS", colnames(dat))]
rna.vars = colnames(dat)[grepl("rna.PC", colnames(dat))]
clin.vars = c("t0_madrs", "sex", "age", "BMI", "t0_diagn_by_age", 
              paste0("t0_cidi_", c("depression_episode", "depression_recurrent", 
                                   "anxiety", "substance")), 
              paste0("t0_pid_", c("negaff", "detach", "psycho", "antago", "disinh")),
              paste0("t0_bsi_", c("soma", "zwan", "unsi", "depr", "angs", "aggr", 
                                  "phob", "para", "psyc")))

## Source nested cross-validation function
source("./Scripts/functions.R")


# 2 ML Analysis Pipeline---------------------------------

# 2.1 Clinical-------------------------------------------

## Define parameters
x = clin.vars

## Run analysis
clinical.output = ml_pipeline(df = dat[ids_for_analysis,], 
                              x = x,
                              y = "bdi_locf_improve",
                              k_outer = 5, 
                              k_inner = 5, 
                              num_repeats = 10, 
                              parallel = TRUE,
                              seed = 10)


# Save results
save(clinical.output, file = "./Results/clinical.output.RData")



# 2.2 Cytokines------------------------------------------

## Define parameters
x = cyto_ref$vars

## Run analysis
cytokine.output = ml_pipeline(df = dat[ids_for_analysis,], 
                              x = x,
                              y = "bdi_locf_improve",
                              k_outer = 5, 
                              k_inner = 5, 
                              num_repeats = 100, 
                              parallel = TRUE,
                              seed = 11)

# Save results
save(cytokine.output, file = "./Results/cytokine.output.RData")




# 2.3 Gene-expression------------------------------------

## Define parameters
x = rna.vars

## Run analysis
rna.output = ml_pipeline(df = dat[ids_for_analysis,], 
                         x = x,
                         y = "bdi_locf_improve",
                         k_outer = 5, 
                         k_inner = 5, 
                         num_repeats = 100, 
                         parallel = TRUE,
                         seed = 12)

# Save results
save(rna.output, file = "./Results/rna.output.RData")




# 2.4 Polygenic Risk Scores------------------------------


## Define parameters
x = prs.vars

## Run analysis
prs.output = ml_pipeline(df = dat[ids_for_analysis,], 
                         x = x,
                         y = "bdi_locf_improve",
                         k_outer = 5, 
                         k_inner = 5, 
                         num_repeats = 100, 
                         parallel = TRUE,
                         seed = 13)

# Save results
save(prs.output, file = "./Results/prs.output.RData")


# 2.5 Combined Immunophenotyping-------------------------


# 2.5.1 Without Covariates-------------------------------


## Define parameters
x = c(prs.vars, cyto_ref$vars, rna.vars)

## Run analysis
omics.output = ml_pipeline(df = dat[ids_for_analysis,], 
                           x = x,
                           y = "bdi_locf_improve",
                           k_outer = 5, 
                           k_inner = 5, 
                           num_repeats = 100, 
                           parallel = TRUE,
                           seed = 14)

# Save results
save(omics.output, file = "./Results/omics.output.RData")



# 2.5.2 With Covariates----------------------------------

## Define parameters
x = c(clin.vars, prs.vars, cyto_ref$vars, rna.vars)

## Run analysis
omicsplusclin.output = ml_pipeline(df = dat[ids_for_analysis,], 
                                   x = x,
                                   y = "bdi_locf_improve",
                                   k_outer = 5, 
                                   k_inner = 5, 
                                   num_repeats = 100, 
                                   parallel = TRUE,
                                   seed = 15)

# Save results
save(omicsplus.output, file = "./Results/omicsplusclin.output.RData")


# -------------------------------------------------------
# Analysis-----------------------------------------------
# -------------------------------------------------------


# 1 Preparation------------------------------------------

## Load packages
library("tidyverse")
library("caret")
library("glmnet")

## Set Slurm directory
setwd("/home/nkappelmann/OPTIMA/ImmunePrediction")


## Load data
load("./Data/OPTIMA_Cytokine_PreprocessedData.RData")

## Index participants used in analysis
# Index participants with minimum of two BDI observations
ids_with_bdi = dat[, paste0("t", 0:7, "_bdi")] %>% is.na() %>% rowSums < 7

# Index participants with genetic data
ids_with_geneticdata = !is.na(dat$BMI_PRS)

# Show overlap
table(ids_with_bdi, ids_with_geneticdata)

# Index participants with transcriptomic data
ids_with_rnadata = !is.na(dat$rna.PC1) # Present in all individuals


# Index participants with all relevant data
ids_for_analysis = ids_with_bdi & ids_with_geneticdata & ids_with_rnadata

## Load cytokine reference
load("./Data/Cytokine_Reference.RData")


## Define prs.vars
prs.vars = colnames(dat)[grepl("PRS", colnames(dat))]


## Define rna.vars
rna.vars = colnames(dat)[grepl("rna.PC", colnames(dat))]


## Source nested cross-validation function
source("./Scripts/functions.R")


# 2 ML Analysis Pipeline---------------------------------

# 2.1 Baseline Covariates--------------------------------

## Define parameters
x = c("t0_bdi_std", "sex_std", "age_std", "BMI_std")

## Run analysis
covariates.output = nested.cv(data = dat[ids_for_analysis,], 
                              x = x,
                              y = "t7_bdi_locf",
                              k.outer = 5, 
                              k.inner = 5, 
                              num_repeats = 100, 
                              runGLMnet = TRUE,
                              runRF = TRUE,
                              runKNN = TRUE,
                              seed = 10)

# Save results
save(covariates.output, file = "./Results/covariates.output.RData")



# 2.2 Cytokines------------------------------------------

## Define parameters
x = cyto_ref$vars

## Run analysis
cytokine.output = nested.cv(data = dat[ids_for_analysis,], 
                            x = x,
                            y = "t7_bdi_locf",
                            k.outer = 5, 
                            k.inner = 5, 
                            num_repeats = 100, 
                            runGLMnet = TRUE,
                            runRF = TRUE,
                            runKNN = TRUE,
                            seed = 11)

# Save results
save(cytokine.output, file = "./Results/cytokine.output.RData")




# 2.3 Gene-expression------------------------------------

## Define parameters
x = rna.vars

## Run analysis
rna.output = nested.cv(data = dat[ids_for_analysis,], 
                       x = x,
                       y = "t7_bdi_locf",
                       k.outer = 5, 
                       k.inner = 5, 
                       num_repeats = 100, 
                       runGLMnet = TRUE,
                       runRF = TRUE,
                       runKNN = TRUE,
                       seed = 12)

# Save results
save(rna.output, file = "./Results/rna.output.RData")




# 2.4 Polygenic Risk Scores------------------------------


## Define parameters
x = prs.vars

## Run analysis
prs.output = nested.cv(data = dat[ids_for_analysis,], 
                       x = x,
                       y = "t7_bdi_locf",
                       k.outer = 5, 
                       k.inner = 5, 
                       num_repeats = 100, 
                       runGLMnet = TRUE,
                       runRF = TRUE,
                       runKNN = TRUE,
                       seed = 13)

# Save results
save(prs.output, file = "./Results/prs.output.RData")


# 2.5 Combined Immunophenotyping-------------------------


# 2.5.1 Without Covariates-------------------------------


## Define parameters
x = c(prs.vars, cyto_ref$vars, rna.vars)

## Run analysis
omics.output = nested.cv(data = dat[ids_for_analysis,], 
                         x = x,
                         y = "t7_bdi_locf",
                         k.outer = 5, 
                         k.inner = 5, 
                         num_repeats = 100, 
                         runGLMnet = TRUE,
                         runRF = TRUE,
                         runKNN = TRUE,
                         seed = 14)

# Save results
save(omics.output, file = "./Results/omics.output.RData")



# 2.5.2 With Covariates----------------------------------

## Define parameters
x = c(prs.vars, cyto_ref$vars, rna.vars, "t0_bdi_std", "sex_std", "age_std", "BMI_std")

## Run analysis
omicsplus.output = nested.cv(data = dat[ids_for_analysis,], 
                            x = x,
                            y = "t7_bdi_locf",
                            k.outer = 5, 
                            k.inner = 5, 
                            num_repeats = 100, 
                            runGLMnet = TRUE,
                            runRF = TRUE,
                            runKNN = TRUE,
                            seed = 15)

# Save results
save(omicsplus.output, file = "./Results/omicsplus.output.RData")


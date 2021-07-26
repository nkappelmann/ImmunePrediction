

# -------------------------------------------------------
# Analysis-----------------------------------------------
# -------------------------------------------------------


# 1 Preparation------------------------------------------

## Load packages
library("tidyverse")
library("ggcorrplot")
library("RColorBrewer")
library("gghalves")
library("caret")
library("glmnet")
library("PubHelper")

## Set Slurm directory
setwd("/home/nkappelmann/OPTIMA/ImmunePrediction")


## Load data
load("./Data/OPTIMA_Cytokine_PreprocessedData.RData")

## Index participants used in analysis
# Index participants with minimum of two BDI observations
ids_with_bdi = dat[, paste0("t", 0:7, "_bdi")] %>% is.na() %>% rowSums < 7

# Index participants with genetic data
ids_with_geneticdata = !is.na(dat$BMI_PRS)

# Show overview
table(ids_with_bdi, ids_with_geneticdata)


# Index participants with transcriptomic data
ids_with_rnadata = !is.na(dat$ENSG00000130066) # Present in all individuals


# Index participants with all relevant data
ids_for_analysis = ids_with_bdi & ids_with_geneticdata & ids_with_rnadata

## Load cytokine reference
load("./Data/Cytokine_Reference.RData")


## Define prs.vars
prs.vars = paste0(list.files("/home/nkappelmann/PRScs/ImmunePrediction/"), 
                  "_PRS")



## Define rna.vars
rna.vars = colnames(dat)[grepl("ENSG", colnames(dat))]


## Source nested cross-validation function
source("./Scripts/functions.R")


# 3 Run analyses-----------------------------------------

# 3.1 Covariates-----------------------------------------

## Define parameters
x = c("t0_bdi_std", "sex_std", "age_std", "BMI_std")

## Run analysis
set.seed(8)
covariates.base.output = nested.cv(data = dat[ids_for_analysis,], 
                                   x = x,
                                   y = "t7_bdi_locf",
                                   k.outer = 5, 
                                   k.inner = 5, 
                                   num_repeats = 100, 
                                   perm.test = FALSE,
                                   runGLMnet = TRUE,
                                   runRF = TRUE,
                                   runKNN = TRUE,
                                   runNNET = FALSE,
                                   runSVM = FALSE,
                                   runBART = FALSE)

# Save results
save(covariates.base.output, file = "./Results/covariates.base.output.RData")

## Run permutation analysis
set.seed(9)
covariates.base.perm = nested.cv(data = dat[ids_for_analysis,], 
                                 x = x,
                                 y = "t7_bdi_locf",
                                 k.outer = 5, 
                                 k.inner = 5, 
                                 num_repeats = 100, 
                                 perm.test = TRUE,
                                 runGLMnet = TRUE,
                                 runRF = TRUE,
                                 runKNN = TRUE,
                                 runNNET = FALSE,
                                 runSVM = FALSE,
                                 runBART = FALSE)

save(covariates.base.perm, file = "./Results/covariates.base.perm.RData")




# 3.2 Cytokines------------------------------------------


# 3.2.1 Without Covariates-------------------------------

## Define parameters
x = cyto_ref$vars

## Run analysis
set.seed(10)
cytokine.base.output = nested.cv(data = dat[ids_for_analysis,], 
                                 x = x,
                                 y = "t7_bdi_locf",
                                 k.outer = 5, 
                                 k.inner = 5, 
                                 num_repeats = 100, 
                                 perm.test = FALSE,
                                 runGLMnet = TRUE,
                                 runRF = TRUE,
                                 runKNN = TRUE,
                                 runNNET = FALSE,
                                 runSVM = FALSE,
                                 runBART = FALSE)

# Save results
save(cytokine.base.output, file = "./Results/cytokine.base.output.RData")

## Run permutation analysis
set.seed(11)
cytokine.base.perm = nested.cv(data = dat[ids_for_analysis,], 
                               x = x,
                               y = "t7_bdi_locf",
                               k.outer = 5, 
                               k.inner = 5, 
                               num_repeats = 100, 
                               perm.test = TRUE,
                               runGLMnet = TRUE,
                               runRF = TRUE,
                               runKNN = TRUE,
                               runNNET = FALSE,
                               runSVM = FALSE,
                               runBART = FALSE)

save(cytokine.base.perm, file = "./Results/cytokine.base.perm.RData")



# 3.2.2 With Covariates----------------------------------


## Define parameters
x = c(cyto_ref$vars, "t0_bdi_std", "sex_std", "age_std", "BMI_std")

## Run analysis
set.seed(12)
cytokine.comb.output = nested.cv(data = dat[ids_for_analysis,], 
                                 x = x,
                                 y = "t7_bdi_locf",
                                 k.outer = 5, 
                                 k.inner = 5, 
                                 num_repeats = 100, 
                                 perm.test = FALSE,
                                 runGLMnet = TRUE,
                                 runRF = TRUE,
                                 runKNN = TRUE,
                                 runNNET = FALSE,
                                 runSVM = FALSE,
                                 runBART = FALSE)

# Save results
save(cytokine.comb.output, file = "./Results/cytokine.comb.output.RData")

## Run permutation analysis
set.seed(13)
cytokine.comb.perm = nested.cv(data = dat[ids_for_analysis,], 
                               x = x,
                               y = "t7_bdi_locf",
                               k.outer = 5, 
                               k.inner = 5, 
                               num_repeats = 100, 
                               perm.test = TRUE,
                               runGLMnet = TRUE,
                               runRF = TRUE,
                               runKNN = TRUE,
                               runNNET = FALSE,
                               runSVM = FALSE,
                               runBART = FALSE)

save(cytokine.comb.perm, file = "./Results/cytokine.comb.perm.RData")



# 3.3 Transcriptomics------------------------------------


# 3.3.1 Without Covariates-------------------------------


## Define parameters
x = rna.vars

## Run analysis
set.seed(10)
rna.base.output = nested.cv(data = dat[ids_for_analysis,], 
                            x = x,
                            y = "t7_bdi_locf",
                            k.outer = 5, 
                            k.inner = 5, 
                            num_repeats = 100, 
                            perm.test = FALSE,
                            runGLMnet = TRUE,
                            runRF = TRUE,
                            runKNN = TRUE,
                            runNNET = FALSE,
                            runSVM = FALSE,
                            runBART = FALSE)

# Save results
save(rna.base.output, file = "./Results/rna.base.output.RData")

## Run permutation analysis
set.seed(11)
rna.base.perm = nested.cv(data = dat[ids_for_analysis,], 
                          x = x,
                          y = "t7_bdi_locf",
                          k.outer = 5, 
                          k.inner = 5, 
                          num_repeats = 100, 
                          perm.test = TRUE,
                          runGLMnet = TRUE,
                          runRF = TRUE,
                          runKNN = TRUE,
                          runNNET = FALSE,
                          runSVM = FALSE,
                          runBART = FALSE)

save(rna.base.perm, file = "./Results/rna.base.perm.RData")


# 3.3.2 With Covariates----------------------------------


## Define parameters
x = c(rna.vars, "t0_bdi_std", "sex_std", "age_std", "BMI_std")

## Run analysis
set.seed(12)
rna.comb.output = nested.cv(data = dat[ids_for_analysis,], 
                            x = x,
                            y = "t7_bdi_locf",
                            k.outer = 5, 
                            k.inner = 5, 
                            num_repeats = 100, 
                            perm.test = FALSE,
                            runGLMnet = TRUE,
                            runRF = TRUE,
                            runKNN = TRUE,
                            runNNET = FALSE,
                            runSVM = FALSE,
                            runBART = FALSE)

# Save results
save(rna.comb.output, file = "./Results/rna.comb.output.RData")

## Run permutation analysis
set.seed(13)
rna.comb.perm = nested.cv(data = dat[ids_for_analysis,], 
                          x = x,
                          y = "t7_bdi_locf",
                          k.outer = 5, 
                          k.inner = 5, 
                          num_repeats = 100, 
                          perm.test = TRUE,
                          runGLMnet = TRUE,
                          runRF = TRUE,
                          runKNN = TRUE,
                          runNNET = FALSE,
                          runSVM = FALSE,
                          runBART = FALSE)

save(rna.comb.perm, file = "./Results/rna.comb.perm.RData")






# 3.4 Polygenic Risk Scores------------------------------


# 3.4.1 Without Covariates-------------------------------


## Define parameters
x = prs.vars

## Run analysis
set.seed(10)
prs.base.output = nested.cv(data = dat[ids_for_analysis,], 
                            x = x,
                            y = "t7_bdi_locf",
                            k.outer = 5, 
                            k.inner = 5, 
                            num_repeats = 100, 
                            perm.test = FALSE,
                            runGLMnet = TRUE,
                            runRF = TRUE,
                            runKNN = TRUE,
                            runNNET = FALSE,
                            runSVM = FALSE,
                            runBART = FALSE)

# Save results
save(prs.base.output, file = "./Results/prs.base.output.RData")

## Run permutation analysis
set.seed(11)
prs.base.perm = nested.cv(data = dat[ids_for_analysis,], 
                          x = x,
                          y = "t7_bdi_locf",
                          k.outer = 5, 
                          k.inner = 5, 
                          num_repeats = 100, 
                          perm.test = TRUE,
                          runGLMnet = TRUE,
                          runRF = TRUE,
                          runKNN = TRUE,
                          runNNET = FALSE,
                          runSVM = FALSE,
                          runBART = FALSE)

save(prs.base.perm, file = "./Results/prs.base.perm.RData")




# 3.2.2 With Covariates----------------------------------



## Define parameters
x = c(prs.vars, "t0_bdi_std", "sex_std", "age_std", "BMI_std")

## Run analysis
set.seed(12)
prs.comb.output = nested.cv(data = dat[ids_for_analysis,], 
                            x = x,
                            y = "t7_bdi_locf",
                            k.outer = 5, 
                            k.inner = 5, 
                            num_repeats = 100, 
                            perm.test = FALSE,
                            runGLMnet = TRUE,
                            runRF = TRUE,
                            runKNN = TRUE,
                            runNNET = FALSE,
                            runSVM = FALSE,
                            runBART = FALSE)

# Save results
save(prs.comb.output, file = "./Results/prs.comb.output.RData")

## Run permutation analysis
set.seed(13)
prs.comb.perm = nested.cv(data = dat[ids_for_analysis,], 
                          x = x,
                          y = "t7_bdi_locf",
                          k.outer = 5, 
                          k.inner = 5, 
                          num_repeats = 100, 
                          perm.test = TRUE,
                          runGLMnet = TRUE,
                          runRF = TRUE,
                          runKNN = TRUE,
                          runNNET = FALSE,
                          runSVM = FALSE,
                          runBART = FALSE)

save(prs.comb.perm, file = "./Results/prs.comb.perm.RData")





# 3.5 Combined Immunophenotyping-------------------------


# 3.5.1 Without Covariates-------------------------------


## Define parameters
x = c(prs.vars, cyto_ref$vars, rna.vars)

## Run analysis
set.seed(10)
omics.base.output = nested.cv(data = dat[ids_for_analysis,], 
                              x = x,
                              y = "t7_bdi_locf",
                              k.outer = 5, 
                              k.inner = 5, 
                              num_repeats = 100, 
                              perm.test = FALSE,
                              runGLMnet = TRUE,
                              runRF = TRUE,
                              runKNN = TRUE,
                              runNNET = FALSE,
                              runSVM = FALSE,
                              runBART = FALSE)

# Save results
save(omics.base.output, file = "./Results/omics.base.output.RData")

## Run permutation analysis
set.seed(11)
omics.base.perm = nested.cv(data = dat[ids_for_analysis,], 
                            x = x,
                            y = "t7_bdi_locf",
                            k.outer = 5, 
                            k.inner = 5, 
                            num_repeats = 100, 
                            perm.test = TRUE,
                            runGLMnet = TRUE,
                            runRF = TRUE,
                            runKNN = TRUE,
                            runNNET = FALSE,
                            runSVM = FALSE,
                            runBART = FALSE)

save(omics.base.perm, file = "./Results/omics.base.perm.RData")



# 3.5.2 With Covariates----------------------------------



## Define parameters
x = c(prs.vars, cyto_ref$vars, rna.vars, "t0_bdi_std", "sex_std", "age_std", "BMI_std")

## Run analysis
set.seed(12)
omics.comb.output = nested.cv(data = dat[ids_for_analysis,], 
                              x = x,
                              y = "t7_bdi_locf",
                              k.outer = 5, 
                              k.inner = 5, 
                              num_repeats = 100, 
                              perm.test = FALSE,
                              runGLMnet = TRUE,
                              runRF = TRUE,
                              runKNN = TRUE,
                              runNNET = FALSE,
                              runSVM = FALSE,
                              runBART = FALSE)

# Save results
save(omics.comb.output, file = "./Results/omics.comb.output.RData")

## Run permutation analysis
set.seed(13)
omics.comb.perm = nested.cv(data = dat[ids_for_analysis,], 
                            x = x,
                            y = "t7_bdi_locf",
                            k.outer = 5, 
                            k.inner = 5, 
                            num_repeats = 100, 
                            perm.test = TRUE,
                            runGLMnet = TRUE,
                            runRF = TRUE,
                            runKNN = TRUE,
                            runNNET = FALSE,
                            runSVM = FALSE,
                            runBART = FALSE)

save(omics.comb.perm, file = "./Results/omics.comb.perm.RData")





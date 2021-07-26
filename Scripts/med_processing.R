# --------------------------------------------
# Medication Coding---------------------------
# --------------------------------------------



# 1 Preparation-------------------------------

## Load packages
library("tidyverse")
library("readxl")


## Load data
med = read_excel("/binder/mgp/datasets/2020_ImmuneDepression/raw_pheno_optima/Medication/Kopie von qry_Medikation_tn1-t1.xlsx")

# Code ID variable
colnames(med)[1] = "ID"

# 2 Coding------------------------------------

# 2.1 Define medication categories------------

## SSRI
ssri = c("cit", "esc", "fluo", "fluv", "paro", "ser")

## SNRI
snri = c("dul", "rebox", "venla")

## Tricyclics (TCA)
tca = c("amioxid", "amitr", "clomi", "dox", "imi", "nortr", "tia", "trimi")

## MAO-inhibitor
mao = c("moclo", "tranyl")

## Atypical antidepressant
atypical = c("ago", "bup", "mirt", "traz")

## Antipsychotic
antipsychotic = c("ari", "ase", "amisul", "benp", "chlor", "cloz", "flup", "halo",
                  "mel", "ola", "pali", "pali", "pera", "perph", "pipa", "prothi", 
                  "quetia", "risp", "zuclo", "zipra")

## Anxiolytic
anxiolytic = c("broma", "dia", "lora", "prega", "oxa", "opi")


## Other
other = c(
      "joh", # Johanneskraut
      "carb", # Carbamazepin
      "lith", # Lithium
      "valp" # Valproic acid
)



# 2.2 Baseline medication---------------------

# SSRI
med$t0_med_ssri = ifelse(rowSums(!is.na(med[, c(paste0("n1_", ssri), 
                                                paste0("0_", ssri))])) > 0,
                         1, 0)

# SNRI
med$t0_med_snri = ifelse(rowSums(!is.na(med[, c(paste0("n1_", snri), 
                                                paste0("0_", snri))])) > 0,
                         1, 0)

# MAO
med$t0_med_mao = ifelse(rowSums(!is.na(med[, c(paste0("n1_", mao), 
                                               paste0("0_", mao))])) > 0,
                        1, 0)


# TCA
med$t0_med_tca = ifelse(rowSums(!is.na(med[, c(paste0("n1_", tca), 
                                               paste0("0_", tca))])) > 0,
                        1, 0)

# Anxiolytic
med$t0_med_anxiolytic = ifelse(rowSums(!is.na(med[, c(paste0("n1_", anxiolytic), 
                                                      paste0("0_", anxiolytic))])) > 0,
                               1, 0)

# Atypical
med$t0_med_atypical = ifelse(rowSums(!is.na(med[, c(paste0("n1_", atypical), 
                                                    paste0("0_", atypical))])) > 0,
                             1, 0)


# Antipsychotic
med$t0_med_antipsychotic = ifelse(rowSums(!is.na(med[, c(paste0("n1_", antipsychotic), 
                                                         paste0("0_", antipsychotic))]))> 0,
                                  1, 0)

# Other
med$t0_med_other = ifelse(rowSums(!is.na(med[, c(paste0("n1_", other), 
                                                 paste0("0_", other))])) > 0,
                          1, 0)


# 3 Save data---------------------------------

save(med, file = "/binder/mgp/datasets/2020_ImmuneDepression/raw_pheno_optima/Medication/Medication_Processed.RData")

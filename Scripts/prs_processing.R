# -------------------------------------
# PRS Processing-----------------------
# -------------------------------------

# 1 Preparation------------------------

## Load packages
library("tidyverse")

## Define directories
prs.dir = "/home/nkappelmann/PRScs/ImmunePrediction/"
out.dir = "/home/nkappelmann/OPTIMA/ImmunePrediction/Data/PRScs/"

## Define phenotypes from folder
pheno = list.files(prs.dir)

## Define function to deal with missing data
completeFun <- function(data, desiredCols) {
      completeVec <- complete.cases(data[, desiredCols])
      return(data[completeVec, ])
}


# 2 Read Files-------------------------

# 2.1 Loop Read-in---------------------

for(i in 1:length(pheno))     {
   
   if(file.exists(paste0(prs.dir, pheno[i], "/", pheno[i], ".sscore"))) {
      ## Read in temporary data
      temp_dat = read.csv(file = paste0(prs.dir, pheno[i], "/", pheno[i], ".sscore"),
                          sep = "", header = TRUE)
      
      # Define prs name
      temp_prs.name = paste0(pheno[i], "_PRS")
      
      # Change colnames
      colnames(temp_dat) = c("FID", "IID", "ALLELE_CT", "NAMED_ALLELE_DOSAGE_SUM",
                             temp_prs.name)
      
      ## Create ID rows for phenotype
      if(i == 1) {prs_dat = data.frame(IID = temp_dat$IID)} 
      
      ## Merge PRS
      prs_dat = merge(prs_dat, temp_dat[, c("IID", temp_prs.name)], 
                      all.x = TRUE, by = "IID")
      
      ## Delete temporary data
      rm(temp_dat)
      rm(temp_prs.name)
   }
}


# 2.2 Save uncorrected PRSs------------

# .csv
write.table(prs_dat, file = paste0(out.dir, "PRScs_uncorrected.csv"), 
            sep = ";", col.names = TRUE, row.names = FALSE, quote = FALSE)

# .RData
save(prs_dat, file = paste0(out.dir, "PRScs_uncorrected.RData"))



# 3 PRS correction---------------------

# 3.1 Load and merge mds components----
mds = read.csv(file = "/binder/common/genotypes/raw_imputations/Become_Optima_Nov2020/01_qc/v1_vs_v2_vs_v3/09_mds2/Become_Optima_no_dups_no_rel_no_outlier_no_het_prunelist_mds.mds",
               header = TRUE, sep = "")

prs_dat = merge(prs_dat, mds[, c("IID", paste0("C", 1:10))], by = "IID", all.x = TRUE)


# 3.2 Residual correction--------------

## Save PRS variables
prs.vars = paste0(pheno, "_PRS")


## TEMP: Print phenotypes that were not available and delete from vector
prs.vars[prs.vars %in% colnames(prs_dat) == FALSE]
prs.vars = prs.vars[prs.vars %in% colnames(prs_dat)]


## Plot mds intercorrelations
png(filename = paste0(out.dir, "mds_corMat.png"))

pairs(prs_dat[, paste0("C", 1:10)], pch = 19)

dev.off()


## Check and save PC intercorrelations
cor_mat = cor(prs_dat[, paste0("C", 1:10)])

write.table(cor_mat, file = paste0(out.dir, "PC_CorMat.csv"), 
            sep = ";", col.names = TRUE, row.names = FALSE, quote = FALSE)


## Residual-based correction based on 10 PCs

# 10 PCs
model = lm(as.matrix(prs_dat[, prs.vars]) ~ prs_dat[, "C1"] + prs_dat[, "C2"] +
                 prs_dat[, "C3"] + prs_dat[, "C4"] +
                 prs_dat[, "C5"] + prs_dat[, "C6"] +
                 prs_dat[, "C7"] + prs_dat[, "C8"] +
                 prs_dat[, "C9"] + prs_dat[, "C10"], 
           data = prs_dat)


## Check correlation of unadjusted and adjusted PRSs
prs_adjCors = diag(cor(prs_dat[, prs.vars], model$residuals)) %>% round(3)

write.table(prs_adjCors, file = paste0(out.dir, "PRS_AdjustCors.csv"), 
            sep = ";", col.names = TRUE, row.names = FALSE, quote = FALSE)



## Overwrite uncorrected variables
prs_dat[, prs.vars] = model$residuals





# 3.3 Save corrected PRSs--------------

# .csv
write.table(prs_dat, file = paste0(out.dir, "PRScs_corrected.csv"), 
            sep = ";", col.names = TRUE, row.names = FALSE, quote = FALSE)

# .RData
save(prs_dat,  file = paste0(out.dir, "PRScs_corrected.RData"))

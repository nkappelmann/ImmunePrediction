# ----------------------------------------
# Processing Transcriptomic Data----------
# ----------------------------------------


# 1 Preparation---------------------------


## Load packages
library("tidyverse")
library("caret")

## Define directory
dir = "/home/nkappelmann/OPTIMA/ImmunePrediction/"

## Load DICE transcriptomic reference data
# Note: Obtained from https://dice-database.org/downloads
dice = read.csv(file = paste0(dir, "Data/mean_tpm_merged.csv"), 
                header = TRUE)


## Load OPTIMA RNAseq data
rna = readRDS(file = "/binder/mgp/workspace/ImmuneDepression/04_RNAseq/03_full_data/02_data/07_postprocessing/rna_filtered_batch_corrected_tpm.rds")

## Load IDs
ids = readRDS(file = "/binder/mgp/workspace/ImmuneDepression/04_RNAseq/03_full_data/02_data/07_postprocessing/matched_ids.Rds")


# 2 DICE Processing-----------------------

# Check #genes
nrow(dice) # 57773

## Create variable to index if gene is expressed in any immune cell
dice$immune_expressed = ifelse(rowSums(dice[, 2:ncol(dice)]) > 0, 1, 0)

# Check #genes expressed
table(dice$immune_expressed) # 49013 genes expressed
prop.table(table(dice$immune_expressed)) * 100 # 85% of genes expressed

# Check #genes expressed in monocytes
table(dice$Monocyte..classical != 0 & dice$Monocyte..non.classical) # 35828
prop.table(table(dice$Monocyte..classical != 0 & dice$Monocyte..non.classical)) * 100 # 62%

# Check #genes in OPTIMA data that are expressed in immune cells
table(dice[dice$immune_expressed == 1, "gene"] %in% row.names(rna)) # 12169
sum(dice[dice$immune_expressed == 1, "gene"] %in% row.names(rna)) / nrow(rna) * 100 # 96%



## Remove genes not expressed in immune cells
dice = dice[dice$immune_expressed == 1,]

# New #genes in DICE
nrow(dice) # 49013

## Calculate maximum gene expression value in any cell type
dice$max_tpm = apply(dice[, 2:(ncol(dice) - 1)], 1, max)

## Reduce data to genes with top 10% expression in any immune cell
dice = arrange(dice, desc(max_tpm))
dice = dice[1:round(nrow(dice) * 0.1),]

# New #genes in DICE
nrow(dice) # 4901


# 3 OPTIMA RNAseq Processing--------------

# 3.1 Subsetting--------------------------

# Check #genes in data
nrow(rna) # 12716

## Look at variance of gene transcription across individuals
rna_sd = apply(rna, 1, sd)

## Extract genes with highest variance (i.e., in top 10%)
rna_genes_highSD = head(sort(rna_sd, decreasing = TRUE), round(0.1 * length(rna_sd)))

# New #genes with highest SD
length(rna_genes_highSD) # 1272

## Subset data to genes to those with:
## 1. Top 10% expression in immune cells
## 2. Top 10% variance in OPTIMA/BeCOME samples
rna = rna[row.names(rna) %in% names(rna_genes_highSD) & row.names(rna) %in% dice$gene,]

# New #genes in data
nrow(rna) # 827


# 3.2 Standardisation---------------------

## Transpose data
rna = as.data.frame(t(rna))

## Scale columns
rna = scale(rna)


# 3.3 PCA---------------------------------

# Note: PCA is done to reduce no. of variables

## Run PCA & save to new object
set.seed(10)
rna.pca = preProcess(rna, method = "pca")

# Extract number of components explaining >=95% of variance
rna.pca$numComp # 87 components

# Save PC data 
rna_proc = predict(rna.pca, newdata = rna)

## Change column names
colnames(rna_proc) = paste0("rna.PC", 1:ncol(rna_proc))


# 4 Merge OPTIMA IDs----------------------

# 4.1 Merge individual gene data----------

## Change format to data.frame
rna = as.data.frame(rna)

## Save RNA-ID and remove "_0" from end of strings
rna$ID_rna = gsub("_0", "", row.names(rna))


## Merge IDs
rna = merge(rna, ids[, c("combined_id", "ID")], by.x = "ID_rna",
            by.y = "combined_id", all.x = TRUE)


# 4.2 Merge PC data-----------------------

## Change format to data.frame
rna_proc = as.data.frame(rna_proc)

## Save RNA-ID and remove "_0" from end of strings
rna_proc$ID_rna = gsub("_0", "", row.names(rna_proc))


## Merge IDs
rna_proc = merge(rna_proc, ids[, c("combined_id", "ID")], by.x = "ID_rna",
            by.y = "combined_id", all.x = TRUE)



# 5 Save Data-----------------------------

save(rna, file = paste0(dir, "Data/TranscriptomicData_processed.RData"))
save(rna_proc, file = paste0(dir, "Data/TranscriptomicData_processed_PConly.RData"))

# ---------------------------------------------------------
# Folkersen SumStats Formatting----------------------------
# ---------------------------------------------------------



# 1 Preparation--------------------------------------------

## Load packages
library("tidyverse")


## Define directories
in.dir = "/home/nkappelmann/GWAS_SummaryStats/Folkersen_raw/"
out.dir = "/home/nkappelmann/GWAS_SummaryStats/Folkersen_processed/"
ref.dir = "/home/nrost/PRS/summary_stats/crp/"


## Define files and phenotypes
gwas.files = list.files(dir)
pheno = gsub(".txt", "", gwas.files)


# 2 Run Formatting loop------------------------------------

for(i in 1:length(gwas.files))    {
      
      ## Read in file
      sumstats = read.csv(file = paste0(in.dir, gwas.files[i]), 
                          header = TRUE, sep = "", stringsAsFactors = FALSE)
      
      
      # 2.1 Preliminary formatting-------------------------
      
      ## Set colnames
      colnames(sumstats) = c("MarkerName", "A1", "A2", "eaf", "FreqSE", 
                             "BETA", "SE", "P", "Direction", "samplesize")
      
      ## Set alleles to uppercase
      sumstats$A1 = toupper(sumstats$A1)
      sumstats$A2 = toupper(sumstats$A2)
      
      ## Delete duplicates
      sumstats = sumstats[!duplicated(sumstats$MarkerName),]
      
      ## Split MarkerName variable
      sumstats[, c("chr", "pos", "alleles")] = sumstats %>% 
            separate(MarkerName, c("chr", "pos", "alleles"), sep = "([:])")
      
      ## Set ID variable
      sumstats$ID = paste(sumstats$chr, sumstats$pos, sep = ":")
      
      
      # 2.2 Append rsid------------------------------------
      
      ## Create rsids column
      sumstats$SNP = NA
      
      ## Loop over chromosomes
      for(j in 1:22)    {
            
            ## Read in chromosome reference data
            chr = read.csv(file = paste0(ref.dir, "SNPs_Chr", j, ".txt"),
                           header = TRUE, sep = "", stringsAsFactors = FALSE)
            
            ## Change column names
            colnames(chr) = c("chr", "start", "end", "SNP")
            
            ## Set chromosome variable to numeric
            chr$chr = j
            
            ## Create ID
            chr$ID = paste(chr$chr, chr$start, sep = ":")
            
            ## Subset reference data to present IDs
            chr = chr[chr$ID %in% sumstats$ID,]
            
            ## Delete duplicates
            chr = chr[!duplicated(chr$ID),]
            
            ## Merge SNP ID
            sumstats = sumstats
            
      }
      
      
      
}

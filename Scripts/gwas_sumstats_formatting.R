# -----------------------------------------------
# GWAS Sumstats Formatting-----------------------
# -----------------------------------------------


# 1 Preparation----------------------------------

## Save input and output directories
in.dir = "/home/nkappelmann/GWAS_SummaryStats/AholaOlli/"
out.dir = "/home/nkappelmann/GWAS_SummaryStats/AholaOlli_PRScs/"

## Save files
files = list.files(in.dir)


# 2 Formatting-----------------------------------


for(i in files)   {
      
      ## Read in data
      temp_dat = read.csv(file = paste0(in.dir, i),
                          header = TRUE,
                          sep = "",
                          stringsAsFactors = FALSE)
      
      ## Change column names
      colnames(temp_dat) = c("SNP", "CHR", "pos", "A2", "A1", 
                             "BETA", "SE", "Direction", "P", "HetPVal")
      
      ## Subset data to relevant columns
      temp_dat = temp_dat[, c("SNP", "A1", "A2", "BETA", "P")]
      
      ## Capitalise alleles
      temp_dat[, c("A1", "A2")] = toupper(temp_dat[, c("A1", "A2")])
      
      ## Exclude ".data" from filename for export
      short_file_name = gsub(".data", "", i)
      
      ## Export data
      write.table(temp_dat, file = paste0(out.dir, short_file_name),
                  quote = FALSE, sep = "\t",
                  row.names = FALSE, col.names = TRUE)
      
}

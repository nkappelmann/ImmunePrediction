# --------------------------------
# Collate PRScs files-------------
# --------------------------------


## Read in files
dir = "/home/nkappelmann/PRScs/ImmunePrediction/"

setwd(dir)

pheno = list.files()

# temporary: pheno = pheno[pheno != "PDFGFbb"]

for(i in pheno)   {
      
      setwd(paste0(dir, i))
      
      temp.files = list.files()
      
      for(j in 1:length(temp.files))    {
            
            if(j == 1)  {
                  txt = read.csv(file = paste0(getwd(), "/", temp.files[j]), 
                                 header = FALSE, stringsAsFactors = FALSE,
                                 sep = "") 
            } else      {
                  txt = rbind.data.frame(txt,
                                         read.csv(file = paste0(getwd(), "/", temp.files[j]), 
                                                  header = FALSE, stringsAsFactors = FALSE,
                                                  sep = ""))
            }
            
      }
      
      write.table(txt, file = paste0(i, ".txt"), col.names = FALSE, 
                  row.names = FALSE, quote = FALSE, sep = "\t")
      
      
}
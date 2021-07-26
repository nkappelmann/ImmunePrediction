
#######
### Function to create folder with different seed value
#######
	


folder <- function(folder.name,seed){
	
	FN <- substitute(folder.name)
	FN <- as.character(FN)
	
	if(FN == "dataset")
	{
		FN <- substitute(folder.name)
		FN <- as.character(FN)
		xy <- paste(getwd(),"/",FN,"_",seed, sep="")
		dir.create(xy)
		setwd(paste(getwd(),"/",FN,"_",seed, sep=""))
		}
		else
		{
						
		FN <- substitute(folder.name)
		FN <- as.character(FN)
		xy <- paste(getwd(),"/", FN,"_",seed, sep="")
		dir.create(xy)
		}
}

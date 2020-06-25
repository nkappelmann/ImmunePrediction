# ------------------------------------------------
# Functions---------------------------------------
# ------------------------------------------------


# 1 Preparation-----------------------------------

## Load Packages
library("tidyverse")


# 2 BaselineTable---------------------------------

baselineTable <- function(data, vars, labels = NULL, grouping.var, round_dec = 2,
                          print.vars = FALSE)     {
      
      # Set labels
      if(is.null(labels)){labels = vars}
      
      # Give error if labels have different length
      if(length(labels) != length(vars))  {
            stop("labels need to have the same length as vars.")
      }
      
      
      # Infer classes and levels of variables
      var.details = data.frame(var = vars,
                               class = NA,
                               levels = NA)
      for(i in vars) {
            var.details[var.details$var == i, "class"] = class(data[, i])
            var.details[var.details$var == i, "levels"] = length(unique(data[!is.na(data[, i]), i]))
      }
      
      # Define table class
      var.details$output = ifelse(var.details$class %in% c("numeric", "integer") & 
                                        var.details$levels > 2, "cont", "cat")
      
      ## Create output table
      output = createOutputTable(data = data, vars = vars, labels = labels, 
                                 var.details = var.details,
                                 grouping.var = grouping.var, round_dec = round_dec)
      
      if(print.vars == FALSE) {output$vars = NULL}
      
      ## Return output
      return(output)
      
}

# 3 createEmptyOutputTable------------------------


createOutputTable <- function(data, vars, labels, var.details, grouping.var, round_dec)   {
      
      ## Create empty table
      output = data.frame(vars = character(), 
                          description = character(),
                          statistic = character())
      
      ## Get nrow
      nrow_data = nrow(data)
      
      for(i in 1:length(vars))    {
            
            ## Get new row numbers
            output.type = var.details[i, "output"]
            newrows_start = nrow(output) + 1
            newrows_end = newrows_start + ifelse(output.type == "cont", 3, 
                                 var.details[i, "levels"] + 1)
            newrows = newrows_start:newrows_end
            
            ## Get unique values if categorical output.type
            if(output.type == "cat")      {
                  output.values = unique(data[, vars[i]])
                  output.values = sort(output.values[!is.na(output.values)])
            }
            
            ## Create new empty rows
            output[newrows,] = NA
            
            ## Set vars and labels
            output[newrows, "vars"] = vars[i]
            output[newrows_start, "description"] = labels[i]
            output[newrows_end, "description"] = "   N Missing (%)"
            
            ## Write differently for continuous and categorical variables
            if(output.type == "cont")     {
                  # Set descriptions
                  output[newrows[2:3], "description"] = paste0("   ", c("Mean (SD)",
                                                                        "Median (IQR)"))
                  
                  # Fill Mean (SD)
                  output[newrows[2], "statistic"] = 
                        paste0(round(mean(data[,vars[i]], na.rm = TRUE), round_dec), " (",
                               round(sd(data[, vars[i]], na.rm = TRUE), round_dec), ")")
                  
                  # Fill Median (IQR)
                  output[newrows[3], "statistic"] = 
                        paste0(round(median(data[,vars[i]], na.rm = TRUE), round_dec), " (",
                               round(summary(data[, vars[i]])[2], round_dec), "-",
                               round(summary(data[, vars[i]])[5], round_dec), ")")
                  
            } else      {
                  
                  # Set descriptions
                  output[newrows[2:(length(newrows) - 1)], "description"] = 
                        as.character(output.values)
                  
                  # Fill N (%)
                  for(j in output.values) {
                        output[output$description == j, "statistic"] = 
                             paste0(sum(data[!is.na(data[, vars[i]]), vars[i]] == j), " (",
                                    round(sum(data[!is.na(data[, vars[i]]), vars[i]] == j) / 
                                                nrow_data * 100, round_dec),
                                    "%)")
                  }
                  
                  # Add zeroes to output.value descriptions
                  output[newrows[2:(length(newrows) - 1)], "description"] = 
                        paste0("   ", output.values)
            }
            
            ## Set missing N (%)
            output[newrows_end, "statistic"] = 
                  paste0(sum(is.na(data[, vars[i]])), " (",
                         round(sum(is.na(data[, vars[i]])) / nrow_data * 100, round_dec), "%)")
            
      }
      
      ## Return output
      return(output)
      
}


# 4 Test function---------------------------------

## Load test data
#data(cars)
#cars$speed_cat = ifelse(cars$speed > 15, "Fast Car", "Slow Car")

# Test
#baselineTable(data = cars, vars = c("dist", "speed", "speed_cat"),
#              labels = c("Stopping distance (ft)", "Speed (mph)", "Speed (group)"))


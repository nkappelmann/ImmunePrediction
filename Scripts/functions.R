# -------------------------------------------------------
# Functions----------------------------------------------
# -------------------------------------------------------


# nested.cv----------------------------------------------

nested.cv <- function(
   data, 
   x,
   y,
   k.outer = 10, # number of outer CV folds
   k.inner = 10, # number of inner CV folds
   num_repeats = 1, # number of CV repeats
   runGLMnet = TRUE,
   runRF = TRUE,
   runKNN = TRUE,
   seed = 1 # set Seed for analyses
)  {
   
   ## Subset data to relevant variables only
   data = data[, c(x, y)]
   
   ## Set running seed, which is incremented by 1 every time it is used
   runningSeed = seed
   
   ## Create folds
   obs.all = nrow(data)
   set.seed(runningSeed)
   fold.pool.outer = rep_len(sample(1:k.outer), length.out = obs.all)
   runningSeed = runningSeed + 1
   
   ## Define fitControl object for caret
   fitControl = trainControl(method = "cv",
                             number = k.inner)
   
   ## Define used models
   used.models = c()
   if(runGLMnet)  {used.models = c(used.models, "glmnet")}
   if(runRF)  {used.models = c(used.models, "rf")}
   if(runKNN)  {used.models = c(used.models, "knn")}
   
   
   ## Define tuneGrids
   glmnet.tuneGrid = expand.grid(alpha = seq(from = 0, to = 1, by = 0.1),
                                 lambda = seq(from = 0, to = 1, by = 0.1))
   rf.tuneGrid = data.frame(mtry = unique(round(seq(from = 2, to = length(x), 
                                                    length.out = 5))))
   knn.tuneGrid = expand.grid(k = 1:25)
   
   
   ## Create output data.frame for fit statistics
   fit.stats = expand.grid(num_repeat = 1:num_repeats,
                           k = 1:k.outer,
                           model = used.models,
                           type = c("pred", "perm"))
   fit.stats[, c("RMSE", "Rsquared", "MAE", "alpha", "lambda", "mtry", "k.knn")] = NA
   
   
   ## Create output data.frame for variable importance
   varImp.stats = fit.stats[, c("num_repeat", "k", "model", "type")]
   varImp.stats[, sort(x)] = NA
   
   
   ## Create output data.frame for individual predictions
   # Save rowID
   data$rowID = paste0("r", 1:nrow(data))
   
   # Create data.frame
   pred.stats = expand.grid(num_repeat = 1:num_repeats,
                            k = 1:k.outer,
                            model = used.models,
                            type = c("pred", "perm"))
   pred.stats[, unique(data$rowID)] = NA
   
   
   
   ## Run repeats
   for(repeats in 1:num_repeats)  {
      
      ## Print repeat number to index progress
      cat(paste0("\nRepeat no.:\t\t", repeats))
      
      ## Y is permuted for permutation test
      set.seed(runningSeed)
      data[, "y_perm"] = sample(data[, y], replace = FALSE)
      runningSeed = runningSeed + 1
      
      ## Assign folds
      set.seed(runningSeed)
      data$fold.outer = sample(fold.pool.outer)
      runningSeed = runningSeed + 1
      
      ## Run outer CV
      for(outer in 1:k.outer) {
         
         ## Print outer fold number to index progress
         cat(paste0("\n---Outer CV-Fold no.:\t", outer))
         
         ## Divide data into train and test sets
         train = data[data$fold.outer != outer,]
         test = data[data$fold.outer == outer,]
         
         # 1 glmnet-------------------------------------------
         
         if(runGLMnet)  {
            
            ## Run inner CV
            set.seed(runningSeed)
            glmnet.fit.pred = train(x = train[,x], y = train[,y], 
                               method = "glmnet", metric = "RMSE", 
                               trControl = fitControl,
                               tuneGrid = glmnet.tuneGrid)
            runningSeed = runningSeed + 1
            
            set.seed(runningSeed)
            glmnet.fit.perm = train(x = train[,x], y = train[,"y_perm"], 
                                    method = "glmnet", metric = "RMSE", 
                                    trControl = fitControl,
                                    tuneGrid = glmnet.tuneGrid)
            runningSeed = runningSeed + 1
            
            ## Save index rows to save fit.stats output
            fit.stats.index.pred = which(fit.stats$num_repeat == repeats & 
                                            fit.stats$k == outer & fit.stats$type == "pred" &
                                            fit.stats$model == "glmnet")
            fit.stats.index.perm = which(fit.stats$num_repeat == repeats & 
                                            fit.stats$k == outer & fit.stats$type == "perm" &
                                            fit.stats$model == "glmnet")
            
            ## Save final tuning parameters
            fit.stats[fit.stats.index.pred, c("alpha", "lambda")] = glmnet.fit.pred$bestTune
            fit.stats[fit.stats.index.perm, c("alpha", "lambda")] = glmnet.fit.perm$bestTune
            
            
            ## Predict in independent test set
            glmnet.preds.pred = predict(glmnet.fit.pred, newdata = test)
            glmnet.preds.perm = predict(glmnet.fit.perm, newdata = test)
            
            # Define pred.stats indices to save prediction results
            pred.stats.index.pred = which(pred.stats$num_repeat == repeats & pred.stats$model == "glmnet" &
                                             pred.stats$k == outer & pred.stats$type == "pred")
            pred.stats.index.perm = which(pred.stats$num_repeat == repeats & pred.stats$model == "glmnet" &
                                             pred.stats$k == outer & pred.stats$type == "perm")
            
            # Save predictions
            pred.stats[pred.stats.index.pred, test$rowID] = glmnet.preds.pred
            pred.stats[pred.stats.index.perm, test$rowID] = glmnet.preds.perm
            
            # Save fit statistics
            fit.stats[fit.stats.index.pred, c("RMSE", "Rsquared", "MAE")] = 
               postResample(test[, y], glmnet.preds.pred)
            fit.stats[fit.stats.index.perm, c("RMSE", "Rsquared", "MAE")] = 
               postResample(test[, y], glmnet.preds.perm)
            
            ## Save variable importance
            glmnet.varImp.pred = varImp(glmnet.fit.pred)$importance
            glmnet.varImp.perm = varImp(glmnet.fit.perm)$importance
            
            # save variable label and order by label
            glmnet.varImp.pred$vars = row.names(glmnet.varImp.pred)
            glmnet.varImp.pred = arrange(glmnet.varImp.pred, vars)
            glmnet.varImp.perm$vars = row.names(glmnet.varImp.perm)
            glmnet.varImp.perm = arrange(glmnet.varImp.perm, vars)
            
            # Save importance values
            varImp.stats[fit.stats.index.pred, sort(x)] = glmnet.varImp.pred$Overall
            varImp.stats[fit.stats.index.perm, sort(x)] = glmnet.varImp.perm$Overall
            
         }
         
         
         # 2 rf-----------------------------------------------
         
         if(runRF == TRUE) {
            
            ## Run inner CV
            set.seed(runningSeed)
            rf.fit.pred = train(x = train[,x], y = train[,y], 
                                method = "rf", metric = "RMSE", 
                                trControl = fitControl,
                                tuneGrid = rf.tuneGrid)
            runningSeed = runningSeed + 1
            
            set.seed(runningSeed)
            rf.fit.perm = train(x = train[,x], y = train[,"y_perm"], 
                                method = "rf", metric = "RMSE", 
                                trControl = fitControl,
                                tuneGrid = rf.tuneGrid)
            runningSeed = runningSeed + 1
            
            ## Save index row to save fit.stats output
            fit.stats.index.pred = which(fit.stats$num_repeat == repeats & 
                                            fit.stats$k == outer & fit.stats$type == "pred" &
                                            fit.stats$model == "rf")
            fit.stats.index.perm = which(fit.stats$num_repeat == repeats & 
                                            fit.stats$k == outer & fit.stats$type == "perm" &
                                            fit.stats$model == "rf")
            
            ## Save final tuning parameters
            fit.stats[fit.stats.index.pred, "mtry"] = rf.fit.pred$finalModel$tuneValue
            fit.stats[fit.stats.index.perm, "mtry"] = rf.fit.perm$finalModel$tuneValue
            
            ## Predict in independent test set
            rf.preds.pred = predict(rf.fit.pred, newdata = test)
            rf.preds.perm = predict(rf.fit.perm, newdata = test)
            
            # Define pred.stats indices to save prediction results
            pred.stats.index.pred = which(pred.stats$num_repeat == repeats & pred.stats$model == "rf" &
                                             pred.stats$k == outer & pred.stats$type == "pred")
            pred.stats.index.perm = which(pred.stats$num_repeat == repeats & pred.stats$model == "rf" &
                                             pred.stats$k == outer & pred.stats$type == "perm")
            
            
            # Save predictions
            pred.stats[pred.stats.index.pred, test$rowID] = rf.preds.pred
            pred.stats[pred.stats.index.perm, test$rowID] = rf.preds.perm
            
            
            # Save fit statistics
            fit.stats[fit.stats.index.pred, c("RMSE", "Rsquared", "MAE")] = 
               postResample(test[, y], rf.preds.pred)
            fit.stats[fit.stats.index.perm, c("RMSE", "Rsquared", "MAE")] = 
               postResample(test[, y], rf.preds.perm)
            
            ## Save variable importance
            rf.varImp.pred = varImp(rf.fit.pred)$importance
            rf.varImp.perm = varImp(rf.fit.perm)$importance
            
            # save variable label and order by label
            rf.varImp.pred$vars = row.names(rf.varImp.pred)
            rf.varImp.pred = arrange(rf.varImp.pred, vars)
            
            rf.varImp.perm$vars = row.names(rf.varImp.perm)
            rf.varImp.perm = arrange(rf.varImp.perm, vars)
            
            # Save importance values
            varImp.stats[fit.stats.index.pred, sort(x)] = rf.varImp.pred$Overall
            varImp.stats[fit.stats.index.perm, sort(x)] = rf.varImp.perm$Overall
            
         }
         
         # 3 knn----------------------------------------------
         
         if(runKNN == TRUE) {
            
            ## Run inner CV
            set.seed(runningSeed)
            knn.fit.pred = train(x = train[,x], y = train[,y], 
                                 method = "knn", metric = "RMSE", 
                                 trControl = fitControl,
                                 tuneGrid = knn.tuneGrid)
            runningSeed = runningSeed + 1
            
            set.seed(runningSeed)
            knn.fit.perm = train(x = train[,x], y = train[,"y_perm"], 
                                 method = "knn", metric = "RMSE", 
                                 trControl = fitControl,
                                 tuneGrid = knn.tuneGrid)
            runningSeed = runningSeed + 1
            
            ## Save index row to save fit.stats output
            fit.stats.index.pred = which(fit.stats$num_repeat == repeats & 
                                            fit.stats$k == outer & fit.stats$type == "pred" &
                                            fit.stats$model == "knn")
            fit.stats.index.perm = which(fit.stats$num_repeat == repeats & 
                                            fit.stats$k == outer & fit.stats$type == "perm" &
                                            fit.stats$model == "knn")
            
            ## Save final tuning parameters
            fit.stats[fit.stats.index.pred, "k.knn"] = knn.fit.pred$finalModel$tuneValue
            fit.stats[fit.stats.index.perm, "k.knn"] = knn.fit.perm$finalModel$tuneValue
            
            ## Predict in independent test set
            knn.preds.pred = predict(knn.fit.pred, newdata = test)
            knn.preds.perm = predict(knn.fit.perm, newdata = test)
            
            # Define pred.stats indices to save prediction results
            pred.stats.index.pred = which(pred.stats$num_repeat == repeats & pred.stats$model == "knn" &
                                             pred.stats$k == outer & pred.stats$type == "pred")
            pred.stats.index.perm = which(pred.stats$num_repeat == repeats & pred.stats$model == "knn" &
                                             pred.stats$k == outer & pred.stats$type == "perm")
            
            # Save predictions
            pred.stats[pred.stats.index.pred, test$rowID] = knn.preds.pred
            pred.stats[pred.stats.index.perm, test$rowID] = knn.preds.perm
            
            
            # Save fit statistics
            fit.stats[fit.stats.index.pred, c("RMSE", "Rsquared", "MAE")] = 
               postResample(test[, y], knn.preds.pred)
            fit.stats[fit.stats.index.perm, c("RMSE", "Rsquared", "MAE")] = 
               postResample(test[, y], knn.preds.perm)
            
            ## Save variable importance
            knn.varImp.pred = varImp(knn.fit.pred)$importance
            knn.varImp.perm = varImp(knn.fit.perm)$importance
            
            # save variable label and order by label
            knn.varImp.pred$vars = row.names(knn.varImp.pred)
            knn.varImp.pred = arrange(knn.varImp.pred, vars)
            
            knn.varImp.perm$vars = row.names(knn.varImp.perm)
            knn.varImp.perm = arrange(knn.varImp.perm, vars)
            
            # Save importance values
            varImp.stats[fit.stats.index.pred, sort(x)] = knn.varImp.pred$Overall
            varImp.stats[fit.stats.index.perm, sort(x)] = knn.varImp.perm$Overall
            
         }
         
         
      }
      
   }
   
   ## Subset fit.stats to fit.stats with available data
   for(i in colnames(fit.stats)[7:ncol(fit.stats)])   {
      if(sum(is.na(fit.stats[, i])) == nrow(fit.stats))  {fit.stats[, i] = NULL}
   }
   

   ## Collate output for fit statistics and variable importance in list
   output = list(fit = fit.stats, 
                 varImp = varImp.stats, 
                 pred = pred.stats)
   
   return(output)
   
}



# getmode------------------------------------------------

getmode <- function(v) {
   
   uniqv = unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
   
}





# aggregate.varImp---------------------------------------


aggregate.varImp <- function(
   obj # nested.cv output object
)   {
   
   ## Define input and output
   input = obj$varImp
   temp = data.frame(num_repeat = numeric(),
                     k = numeric(),
                     model = character(),
                     var = character(),
                     varImp = numeric())
   
   
   ## Fill varImp data.frame
   for(i in 4:ncol(input))   {
      for(j in 1:nrow(input))   {
         
         ## Create empty new row
         temp[(nrow(temp) + 1), ] = NA
         
         ## Fill general info
         temp[nrow(temp), c("num_repeat", "k")] = 
            input[j, c("num_repeat", "k")]
         temp[nrow(temp), "model"] = 
            as.character(input[j, "model"])
         
         ## Fill var and varImp info
         temp[nrow(temp), "var"] = 
            colnames(input)[i]
         temp[nrow(temp), "varImp"] = 
            input[j, i]
      }
      
   }
   
   ## Reduce data and get median info
   output = rank.aggregate.varImp(x = temp)
   
   
   return(output)
}


# rank.aggregate.varImp----------------------------------

rank.aggregate.varImp = function(
   x # input should be varImp data in long format as provided by aggregate.varImp()
)  {
   
   # Create ranking on variable importance
   varImp.ranking = with(x, by(varImp, var, median, na.rm = TRUE))
   varImp.ranking = data.frame(var = names(varImp.ranking),
                               varImp = as.vector(varImp.ranking)) %>%
      arrange(desc(varImp))
   
   
   # Create levels of individual data
   x$var = factor(x$var, levels = varImp.ranking$var)
   
   # Arrange data by varImp
   x = arrange(x, var)
   
   # Create output list
   output = list(ind = x,
                 ranking = varImp.ranking)
   
   return(output)
}



# extract.fitDiff----------------------------------------

extract.fitDiff = function(
   x = NULL # fit statistics output from nested.xv 
)  {
   
   if(is.null(x)) {stop("X needs to be specified")}
   
   
   output = x$fit
   output$RMSE_delta = NA
   output$R2_delta = NA
   
   ## Loop over CV folds, repetitions and models
   for(i in unique(output$num_repeat))   {
      for(j in unique(output$k))   {
         for(k in unique(output$model))  {
            
            # Get row indices
            index.pred = which(output$num_repeat == i & output$k == j & 
                                  output$model == k & output$type == "pred")
            index.perm = which(output$num_repeat == i & output$k == j & 
                                  output$model == k & output$type == "perm")
            
            # Calculate Delta variables
            output[index.pred, "RMSE_delta"] = output[index.perm, "RMSE"] - output[index.pred, "RMSE"]
            output[index.pred, "R2_delta"] = output[index.pred, "Rsquared"] - output[index.perm, "Rsquared"]
            
         }
      }
   }
   
   return(output[, c("RMSE_delta", "R2_delta")])
   
}



# t.test_loop--------------------------------------------

t.test_loop = function(
   x, # fit statistics dataset
   round_dec = 3
)  {
   
   ## Create empty output data.frame
   output = expand.grid(feature.set = unique(x$feature.set),
                        model = unique(x$model))
   output[, c("RMSE_mean", "RMSE_sd", "R2_median", "R2_25quantile", "R2_75quantile",
              "RMSE_delta_mean", "RMSE_delta_sd", "RMSE_delta_median", "RMSE_delta_5quantile", 
              "RMSE_delta_95quantile", "RMSE", "R2", "RMSE_delta", "tval", "pval")] = NA
   
   ## Loop over combinations
   for(i in unique(output$feature.set))   {
      for(j in unique(output$model))   {
         
         ## Subset data
         x.pred = filter(x, feature.set == i & model == j & type == "pred")
         x.perm = filter(x, feature.set == i & model == j & type == "perm")
         
         ## Conduct paired t-test
         t.test.fit = t.test(x.pred$RMSE, x.perm$RMSE, paired = TRUE)
         
         ## Write output data
         output.index = which(output$feature.set == i & output$model == j)
         
         output[output.index, "RMSE_mean"] = mean(x.pred$RMSE)
         output[output.index, "RMSE_sd"] = sd(x.pred$RMSE)
         output[output.index, "R2_median"] = median(x.pred$Rsquared, na.rm = TRUE)
         output[output.index, "R2_25quantile"] = quantile(x.pred$Rsquared, probs = 0.25, na.rm = TRUE)
         output[output.index, "R2_75quantile"] = quantile(x.pred$Rsquared, probs = 0.75, na.rm = TRUE)
         
         RMSE_delta = x.pred$RMSE - x.perm$RMSE
         output[output.index, "RMSE_delta_mean"] = mean(RMSE_delta)
         output[output.index, "RMSE_delta_sd"] = sd(RMSE_delta)
         output[output.index, "RMSE_delta_median"] = median(RMSE_delta)
         output[output.index, c("RMSE_delta_5quantile", "RMSE_delta_95quantile")] = 
            quantile(RMSE_delta, probs = c(0.05, 0.95))
         
         output[output.index, "tval"] = t.test.fit$statistic
         output[output.index, "pval"] = t.test.fit$p.value
         
      }
   }
   
   
   output[, "RMSE"] = with(output, paste0(round(RMSE_mean, round_dec), " (", 
                                          round(RMSE_sd, round_dec), ")"))
   output[, "R2"] = with(output, paste0(round(R2_median, round_dec), " (", 
                                        round(R2_25quantile, round_dec), "-", 
                                        round(R2_25quantile, round_dec), ")"))
   
   output[, "RMSE_delta"] = with(output, paste0(round(RMSE_delta_mean, round_dec), " (", 
                                          round(RMSE_delta_sd, round_dec), ")"))
   
   return(output)
   
}



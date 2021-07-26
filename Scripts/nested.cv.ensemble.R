# nested.cv.ensemble-------------------------------------

nested.cv.ensemble <- function(
      data, 
      x,
      y,
      k.outer = 10, # number of outer CV folds
      k.inner = 10, # number of inner CV folds
      num_repeats = 1, # number of CV repeats
      perm.test = FALSE, # Should y be shuffled to obtain a permutation distribution?
      runGLMnet = TRUE,
      runRF = TRUE,
      runKNN = TRUE
)  {
      
      ## Subset data to relevant variables only
      data = data[, c(x, y)]
      
      ## Create folds
      obs.all = nrow(data)
      fold.pool.outer = rep_len(sample(1:k.outer), length.out = obs.all)
      
      ## Define used models
      used.models = c()
      if(runGLMnet)  {used.models = c(used.models, "glmnet")}
      if(runRF)  {used.models = c(used.models, "rf")}
      if(runKNN)  {used.models = c(used.models, "knn")}
      
      ## Define tuneGrids
      glmnet.tuneGrid = expand.grid(alpha = seq(from = 0, to = 1, by = 0.2),
                                    lambda = seq(from = 0, to = 1, by = 0.2))
      rf.tuneGrid = data.frame(mtry = unique(round(seq(from = 2, to = length(x), 
                                                       length.out = 5))))
      knn.tuneGrid = expand.grid(k = 1:25)
      
      ## Create output data.frame for fit statistics
      fit.stats = expand.grid(num_repeat = 1:num_repeats,
                              k = 1:k.outer,
                              model = used.models)
      fit.stats[, c("RMSE", "Rsquared", "MAE", "alpha", "lambda", "mtry", "k.knn", 
                    "degree", "scale", "C", "size", "decay",
                    "num_trees", "k.bart", "beta", "nu")] = NA
      
      
      ## Create output data.frame for variable importance
      varImp.stats = fit.stats[, c("num_repeat", "k", "model")]
      varImp.stats[, sort(x)] = NA
      
      
      ## Create output data.frame for individual predictions
      # Save rowID
      data$rowID = paste0("r", 1:nrow(data))
      
      # Create data.frame
      pred.stats = expand.grid(num_repeat = 1:num_repeats,
                               model = used.models)
      pred.stats[, unique(data$rowID)] = NA
      
      
      
      ## Run repeats
      for(repeats in 1:num_repeats)  {
            
            ## Print repeat number to index progress
            cat(paste0("\nRepeat no.:\t\t", repeats))
            
            ## Permute Y if a permutation test is indicated
            if(perm.test == TRUE)   {
                  data[, y] = sample(data[, y], replace = FALSE)
            }
            
            ## Assign folds
            data$fold.outer = sample(fold.pool.outer)
            
            ## Run outer CV
            for(outer in 1:k.outer) {
                  
                  ## Print outer fold number to index progress
                  cat(paste0("\n---Outer CV-Fold no.:\t", outer))
                  
                  ## Divide data into train and test sets
                  train = data[data$fold.outer != outer,]
                  test = data[data$fold.outer == outer,]
                  
                  ## Define fitControl object for caret
                  fitControl = trainControl(method = "cv",
                                            number = k.inner,
                                            savePredictions = "final",
                                            index = createResample(train[,y], 
                                                                   k.inner))
                  
                  
                  # 1 Run individual algorithms in inner CV------------
                  
                  ## glmnet
                  glmnet.fit = train(x = train[,x], y = train[,y], 
                                     method = "glmnet", metric = "RMSE", 
                                     trControl = fitControl,
                                     tuneGrid = glmnet.tuneGrid)
                  
                  
                  ## random forest
                  rf.fit = train(x = train[,x], y = train[,y], 
                                 method = "rf", metric = "RMSE", 
                                 trControl = fitControl,
                                 tuneGrid = rf.tuneGrid)
                  
                  ## k-nearest neighbours
                  knn.fit = train(x = train[,x], y = train[,y], 
                                  method = "knn", metric = "RMSE", 
                                  trControl = fitControl,
                                  tuneGrid = knn.tuneGrid)
                  
                  
                  
                  # 2 Run ensemble model in inner CV-------------------
                  
                  # Note: Some code adapted from supplementary material of Kew & Mitchell (2015)
                  
                  # Create a list of all models
                  message("list of all models")
                  model_list = list(glmnet.fit, rf.fit, knn.fit) 
                  names(model_list) <- sapply(model_list, function(x) x$method)
                  class(model_list) = "caretList"
                  sort(sapply(model_list, function(x) min(x$results$RMSE)))
                  
                  #Make a greedy ensemble - currently can only use RMSE
                  message("make a greedy ensemble")
                  greedy <- caretEnsemble(model_list, iter = 1000L)
                  sort(greedy$weights, decreasing=TRUE)
                  greedy$error
                  
                  tuneList = list(
                        glmnet = caretModelSpec(method = "glmnet",
                                                #tuneLength = 5
                                                tuneGrid = glmnet.tuneGrid
                        ), 
                        rf = caretModelSpec(method = "rf",
                                            #tuneLength = 5
                                            tuneGrid = rf.tuneGrid
                        ), 
                        knn = caretModelSpec(method = "knn",
                                             #tuneLength = 5
                                             tuneGrid = knn.tuneGrid
                        ))
                  
                  model_list = caretList(
                        x = train[,x], y = train[,y], 
                        trControl = fitControl,
                        methodList = c("glmnet", "rf", "knn")#, 
                        #tuneList = tuneList
                  )
                  
                  
                  ## Use greedy algorithm 
                  greedy_ensemble = caretEnsemble(
                        model_list, 
                        metric = "RMSE",
                        trControl = fitControl)
                  summary(greedy_ensemble)
                  
                  
                  ## Extract predictions
                  p <- as.data.frame(predict(model_list, newdata=head(testing)))
                  model_preds = lapply(model_list, predict, newdata = test)
                  ensemble_preds = predict(greedy_ensemble, newdata = test)
                  xyplot(resamples(model_list))
                  modelCor(resamples(model_list))
                  
                  
                  ## Save index row to save fit.stats output
                  fit.stats.index = which(fit.stats$num_repeat == repeats & 
                                                fit.stats$k == outer & 
                                                fit.stats$model == "glmnet")
                  
                  ## Save final tuning parameters
                  fit.stats[fit.stats.index, c("alpha", "lambda")] = glmnet.fit$bestTune
                  
                  
                  
                  ## Predict in independent test set
                  glmnet.preds = predict(glmnet.fit, newdata = test)
                  
                  # Save predictions
                  pred.stats[pred.stats$num_repeat == repeats & pred.stats$model == "glmnet", 
                             test$rowID] = glmnet.preds
                  
                  # Save fit statistics
                  fit.stats[fit.stats.index, c("RMSE", "Rsquared", "MAE")] = 
                        postResample(test[, y], glmnet.preds)
                  
                  ## Save variable importance
                  glmnet.varImp = varImp(glmnet.fit)$importance
                  
                  # save variable label and order by label
                  glmnet.varImp$vars = row.names(glmnet.varImp)
                  glmnet.varImp = arrange(glmnet.varImp, vars)
                  
                  # Save importance values
                  varImp.stats[fit.stats.index, sort(x)] = glmnet.varImp$Overall
                  
                  
                  
                  # 2 rf-----------------------------------------------
                  
                  ## Run inner CV
                  rf.fit = train(x = train[,x], y = train[,y], 
                                 method = "rf", metric = "RMSE", 
                                 trControl = fitControl,
                                 tuneGrid = rf.tuneGrid
                  )
                  
                  ## Save index row to save fit.stats output
                  fit.stats.index = which(fit.stats$num_repeat == repeats & 
                                                fit.stats$k == outer & 
                                                fit.stats$model == "rf")
                  
                  ## Save final tuning parameters
                  fit.stats[fit.stats.index, "mtry"] = rf.fit$finalModel$tuneValue
                  
                  ## Predict in independent test set
                  rf.preds = predict(rf.fit, newdata = test)
                  
                  # Save predictions
                  pred.stats[pred.stats$num_repeat == repeats & pred.stats$model == "rf", 
                             test$rowID] = rf.preds
                  
                  
                  # Save fit statistics
                  fit.stats[fit.stats.index, c("RMSE", "Rsquared", "MAE")] = 
                        postResample(test[, y], rf.preds)
                  
                  ## Save variable importance
                  rf.varImp = varImp(rf.fit)$importance
                  
                  # save variable label and order by label
                  rf.varImp$vars = row.names(rf.varImp)
                  rf.varImp = arrange(rf.varImp, vars)
                  
                  # Save importance values
                  varImp.stats[fit.stats.index, sort(x)] = rf.varImp$Overall
                  
                  
                  # 3 knn----------------------------------------------
                  
                  ## Run inner CV
                  knn.fit = train(x = train[,x], y = train[,y], 
                                  method = "knn", metric = "RMSE", 
                                  trControl = fitControl,
                                  tuneGrid = knn.tuneGrid
                  )
                  
                  ## Save index row to save fit.stats output
                  fit.stats.index = which(fit.stats$num_repeat == repeats & 
                                                fit.stats$k == outer & 
                                                fit.stats$model == "knn")
                  
                  ## Save final tuning parameters
                  fit.stats[fit.stats.index, "k.knn"] = knn.fit$finalModel$tuneValue
                  
                  ## Predict in independent test set
                  knn.preds = predict(knn.fit, newdata = test)
                  
                  # Save predictions
                  pred.stats[pred.stats$num_repeat == repeats & pred.stats$model == "knn", 
                             test$rowID] = knn.preds
                  
                  
                  # Save fit statistics
                  fit.stats[fit.stats.index, c("RMSE", "Rsquared", "MAE")] = 
                        postResample(test[, y], knn.preds)
                  
                  ## Save variable importance
                  knn.varImp = varImp(knn.fit)$importance
                  
                  # save variable label and order by label
                  knn.varImp$vars = row.names(knn.varImp)
                  knn.varImp = arrange(knn.varImp, vars)
                  
                  # Save importance values
                  varImp.stats[fit.stats.index, sort(x)] = knn.varImp$Overall
                  
                  
                  
                  
            }
            
      }
      
      ## Subset fit.stats to fit.stats with available data
      for(i in colnames(fit.stats)[7:ncol(fit.stats)])   {
            if(sum(is.na(fit.stats[, i])) == nrow(fit.stats))  {fit.stats[, i] = NULL}
      }
      
      ## Aggregate pred.stats across models using a voting system
      pred.aggregate.stats = pred.aggregate(preds = pred.stats)
      
      
      ## Collate output for fit statistics and variable importance in list
      output = list(fit = fit.stats, 
                    varImp = varImp.stats, 
                    pred = pred.stats,
                    pred.aggregate = pred.aggregate.stats)
      
      return(output)
      
}



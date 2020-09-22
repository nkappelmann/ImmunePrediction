# -------------------------------------------------
# Reproducible caret prediction problem------------
# -------------------------------------------------


# 1 Setup-----------------------------------------

## Load packages
library("caret")


# 2 Create data-----------------------------------

## Create independent vars
set.seed(5)
dat = data.frame (x1 = rnorm(100, mean = 0, sd = 1),
                  x2 = rnorm(100, mean = 0, sd = 1))

## Create y
set.seed(6)
dat$y = 10 + 0.5*dat$x1 + 0.2*dat$x2 + rnorm(100, mean = 0, sd = 1.5)

## Plot correlations
cor(dat)


## Split into train and test
set.seed(7)
inTrain = sample(1:nrow(dat), size = round((nrow(dat) / 5) * 4), replace = FALSE)

train = dat[inTrain,]
test = dat[-inTrain,]

## Specify variable names
x = c("x1", "x2")
y = "y"

# 3 Build model-----------------------------------

## Define model parameters
fitControl = trainControl(method = "cv",
                          number = 5,
                          savePredictions = FALSE)

## Define tune grid
glmnet.tuneGrid = expand.grid(alpha = seq(from = 0, to = 1, by = 0.2),
                              lambda = seq(from = 0, to = 1, by = 0.2))

## Run inner CV
glmnet.fit = train(x = train[,x], y = train[,y], 
                   method = "glmnet", metric = "RMSE", 
                   trControl = fitControl,
                   tuneGrid = glmnet.tuneGrid)

## Save final tuning parameters
fit.stats = data.frame(alpha = glmnet.fit$finalModel$tuneValue[1,1], 
                       lambda = glmnet.fit$finalModel$tuneValue[1,2]) 


## Predict in independent test set
glmnet.preds = predict(glmnet.fit, newdata = test)


glmnet.coefs = coef(glmnet.fit$finalModel, s = glmnet.fit$bestTune$lambda)
manual.preds = glmnet.coefs[1,] + glmnet.coefs[2,]*test$x1 + glmnet.coefs[3,]*test$x1


# Save fit statistics
fit.stats[1, c("RMSE.predict", "Rsquared.predict", "MAE.predict")] = 
      postResample(test[, y], glmnet.preds)
fit.stats[1, c("RMSE.manual", "Rsquared.manual", "MAE.manual")] = 
      postResample(test[, y], manual.preds)


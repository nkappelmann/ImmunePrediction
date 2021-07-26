#Libraries
library(caret)
library(lattice)
library(reshape)
library(plyr)
library(dismo)
library(raster)
library(sp)
library(rJava)
library(randomForest)
library(nnet)
library(devtools)
#install_github('caretEnsemble', 'zachmayer') #Install zach's caretEnsemble package
library(caretEnsemble)



#Load and split data
data <- input     #   <- read.table(paste(data1,".csv",sep=""),sep=",",header=T)
X <- data[,-(ncol(data))]
Y <- data[,(ncol(data))]

num<-row.names(X)


## seed values
a<-c(102,912,854,1257,459) #from Neetikas random number seeds

## Start for loop for different seed numbers

for(r in 1:length(a)) 
{
set.seed(a[r])  #set seed number


## creating a folder for new seed number

folder(dataset_reg,a[r])

FN <- substitute(dataset_reg)
FN <- as.character(FN)
xy <- paste(getwd(),"/", FN,"_",a[r], sep="")


#starting cross fold validation for loop
for(k in 1:10) 
{


	message("Current CV Fold is:")
	message(k)


	##
	#Train/test split for 10 fold CV
	##
	
	number<-kfold(X, k=10)
	num<-cbind(num,number)
	
	#splits the input data into a training set and external validation set
	#split based on fold number - i.e. each fold uses a different split.
	preds<-NULL
	XX<-NULL
	XX<-which(number==k)
	trainset<-X[-XX,]
	testset<-X[XX,]
	trainMedv <- Y[-XX]
	testMedv <- Y[XX]

	message("Set Up CV folds")
	#Setup CV Folds
	#returnData=FALSE saves some space
	myControl <- trainControl(method='cv',10,1,allowParallel=TRUE,
							  returnResamp='none',returnData=FALSE, 
							  savePredictions=TRUE, verboseIter=TRUE,
							  index=createMultiFolds(trainMedv,10,1))
	PP <- c('pca')

	message("train some models")
	#Train some models
	message("rvm")
	#model1a <- train(trainset, trainMedv, method='rvmLinear', metric = "RMSE", trControl=myControl, tuneLength=20, preProcess=PP)
	model1b <- train(trainset, trainMedv, method='rvmRadial', metric = "RMSE", trControl=myControl, tuneLength=10, preProcess=PP)
	model1c <- train(trainset, trainMedv, method='rvmPoly', metric = "RMSE", trControl=myControl, tuneLength=5, preProcess=PP) #larger tune length slows it down a lot - also makes RMSE worse - overfitting?
	message("SVM")
	#model2a <- train(trainset, trainMedv, method='svmLinear', metric = "RMSE", trControl=myControl, tuneLength=20, preProcess=PP)
	model2b <- train(trainset, trainMedv, method='svmRadial', metric = "RMSE", trControl=myControl, tuneLength=10, preProcess=PP)
	model2c <- train(trainset, trainMedv, method='svmPoly', metric = "RMSE", trControl=myControl, tuneLength=5, preProcess=PP)
	message("Gaussian")
	#model3a <- train(trainset, trainMedv, method='gaussprLinear', metric = "RMSE", trControl=myControl, tuneLength=20, preProcess=PP)
	model3b <- train(trainset, trainMedv, method='gaussprRadial', metric = "RMSE", trControl=myControl, tuneLength=10, preProcess=PP)
	model3c <- train(trainset, trainMedv, method='gaussprPoly', metric = "RMSE", trControl=myControl, tuneLength=5, preProcess=PP)
	message("Neural Networks")
	model4 <- train(trainset, trainMedv, method='nnet', metric = "RMSE", trControl=myControl, tuneGrid=expand.grid(.size=c(1,5,10),.decay=c(0,0.001,0.1)), MaxNWts=10000,linout=TRUE,trace=FALSE, preProcess=PP)
	message("nearest neighbours")
	model5 <- train(trainset, trainMedv, method='knn', metric = "RMSE", trControl=myControl, tuneLength=20, preProcess=PP)
	message("treebag")
	model6 <- train(trainset, trainMedv, method='treebag', metric = "RMSE", trControl=myControl, tuneLength=20, preProcess=PP)
	message("rf")
	model7 <- train(trainset, trainMedv, method='rf', metric = "RMSE", trControl=myControl, tuneLength=10,ntree=1000,importance=TRUE, preProcess=PP)
	message("cubist")
	model8 <- train(trainset, trainMedv, method='cubist', metric = "RMSE", trControl=myControl, tuneLength=30, preProcess=PP)
	message("linear methods")
	model9 <- train(trainset, trainMedv, method='pls', metric = "RMSE", trControl=myControl, tuneLength=20, preProcess=PP)
	model10 <- train(trainset, trainMedv, method='pcr', metric = "RMSE", trControl=myControl, tuneLength=20, preProcess=PP)
	

	#Make a list of all the models
	message("list of all models")
	all.models <- list(model1b,model1c,model2b,model2c,model3b,model3c,model4,model5,model6,model7,model8,model9,model10) #fullset
	names(all.models) <- sapply(all.models, function(x) x$method)
	sort(sapply(all.models, function(x) min(x$results$RMSE)))

	#Make a greedy ensemble - currently can only use RMSE
	message("make a greedy ensemble")
	greedy <- caretEnsemble(all.models, iter=1000L)
	sort(greedy$weights, decreasing=TRUE)
	greedy$error

	#Make a linear regression ensemble
	message("make a linear ensemble")
	linear <- caretStack(all.models, method='glm', trControl=trainControl(method='cv'))
	summary(linear$ens_model$finalModel)
	linear$error

	#output results of training
	message("Outputting Result of training")
	sink(paste(xy,"/output_fold_",k,".txt",sep=""))
	print("Results")
	print("Greedy")
	print(greedy)
	print("Linear")
	print(linear)
	sink()

	sink(paste(xy,"/error_fold_",k,".csv",sep=""))
	print("Greedy Error")
	print(greedy$error)
	print("Greedy Weights")
	print(greedy$weights)
	print("Linear Error")
	print(linear$error)
	sink()

	##1-Predict for external test set:
	message("Predicting for external test set")
	preds <- data.frame(sapply(all.models, predict, newdata=testset))
	preds$ens_greedy <- predict(greedy, newdata=testset)
	preds$ens_linear <- predict(linear, newdata=testset)
	preds$actualvalue <- paste(testMedv)
	
	#1-output predictions of testing
	message("output predictions")
	sink(paste(xy,"/predicitions_fold_",k,".csv",sep=""))
	print(preds)
	sink()
	
}
	
}


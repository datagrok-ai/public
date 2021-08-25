#name: RCaretNBTrain
#meta.mlname: RCaretNB
#meta.mlrole: train
#description: Custom R Train func for Naive Bayes using Caret
#language: r
#input: dataframe df
#input: string predict_column
#output: blob model

#>>> Initialize local environment and install/require necessary packages

library(caret)
library(e1071)
renv::init()
install.packages("naivebayes")
library(naivebayes)
install.packages("MLmetrics")
library(MLmetrics)
#install.packages("klaR") 
#library(klaR)

#>>> Extract/Prepare train features and target variable
trainY <- df[[predict_column]]
trainX <- df[, -grep(predict_column, colnames(df))]

#>>> Prepare parameters for model training
tuneGrid <- expand.grid(laplace = 0,  usekernel = FALSE, adjust=1)
#tuneGrid <- expand.grid(usekernel = TRUE,  fL = 1, adjust=3)  #<-- Uncomment for method='nb'
ctrl <- trainControl(method="repeatedcv", repeats=5, summaryFunction=multiClassSummary,
  classProbs=TRUE, savePredictions=FALSE, returnResamp="final", returnData=FALSE, trim=TRUE)

#>>> Train model
trained_model <- caret::train(
  x=trainX, 
  y=trainY, 
  method='naive_bayes',
#  method='nb', 
  preProcess=c("center", "scale"), 
  trControl=ctrl, 
  tuneGrid = tuneGrid
  )

#>>> Save trained model 
saveRDS(trained_model, file=model)

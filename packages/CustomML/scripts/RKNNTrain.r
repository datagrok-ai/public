#name: RKNNTrain
#meta.mlname: RKNN
#meta.mlrole: train
#description: Custom R Train func for KNN using Caret
#language: r
#input: dataframe df
#input: string predict_column
#input: int n_neighbors {category: Parameters}
#input: int distance = 1 {category: Parameters; range:1-2} [Parameter of Minkowski distance] 
#input: string kernel = rectangular {category: Parameters; choices: ["rectangular", "triangular", "epanechnikov", "biweight", "triweight", "cos", "inv", "gaussian", "optimal"]}  
#output: blob model

#>>> Initialize local environment and install-require necessary packages

library(caret)
library(e1071)
renv::init()
install.packages('kknn')
library(kknn)
install.packages("MLmetrics")
library(MLmetrics)

#>>> Extract/Prepare train features and target variable
trainY <- df[[predict_column]]
trainX <- df[, -grep(predict_column, colnames(df))]
 
#>>> Prepare parameters for model training
tuneGrid <- expand.grid(kmax = n_neighbors,  distance = distance, kernel = kernel) #Parms for KKNN
#tuneGrid <-  expand.grid(k = n_neighbors) <-- Uncomment for KNN
ctrl <- trainControl(method="cv", summaryFunction=multiClassSummary, classProbs=TRUE, 
        savePredictions=FALSE, returnResamp="final", returnData=FALSE, trim=TRUE)

#>>> Train model
trained_model <- caret::train(
  x=trainX, 
  y=trainY, 
  method='kknn', 
#  method='knn',
  preProcess=c("center", "scale"), 
  metric="ROC", 
  trControl=ctrl, 
  tuneGrid = tuneGrid
  )

#>>> Save trained model
saveRDS(trained_model, file=model)

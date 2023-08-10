#name: Predictive Model SVM
#description: Predictive model based on Support vector machine (SVM)
#reference: https://en.wikipedia.org/wiki/Support_vector_machine
#language: r
#tags: demo
#sample: demog.csv
#input: dataframe datasetPredict [Input data set table]
#input: column predict {type:categorical} [Name of column contains values to predict]
#input: dataframe dataset {columns:numerical} [Input data set table]
#input: column_list columns {type:numerical} [Input data set table columns]
#input: bool fillMissing = TRUE [Fill missing values using kNN algorithm]
#input: double perc = 0.5 [Percentage of data that goes to training]
#output: graphics [ROC plot of predicted test part of dataset]

# Libraries
require(caret)
require(ggplot2)
require(plotROC)
require(VIM)
require(RANN)

# Fill missing values
if (fillMissing == TRUE) {
  dataset = kNN(dataset, k=5)[, colnames(dataset)]
}

# Load the data and construct indices to divide it into training and test data sets.
dataset <- data.frame(dataset, datasetPredict[[predict]])
colnames(dataset)[length(dataset)] = predict
trainYCol = predict
trainIndex <- createDataPartition(dataset[[trainYCol]], p=perc, list=FALSE)
trainData <- dataset[trainIndex,]
trainY = trainData[[trainYCol]]
trainX <- trainData[, -grep(trainYCol, colnames(dataset))]

testData <- dataset[-trainIndex,]
testY = testData[[trainYCol]]
testX <- testData[, -grep(trainYCol, colnames(dataset))]

# Setup for cross validation
set.seed(73)
ctrl <- trainControl(method="repeatedcv", repeats=5, summaryFunction=multiClassSummary,
    classProbs=TRUE, savePredictions=TRUE)

# Train SVM model
set.seed(73)
svm.model <- caret::train(x=trainX, y=trainY, method="svmLinear",
    preProcess=c("center","scale"), metric="ROC", trControl=ctrl)

# Apply model
set.seed(73)
predictY = predict(svm.model, testX)
probY = predict(svm.model, testX, type="prob")

# ROC plot
df.plot <- data.frame(D=as.numeric(testY) - 1, M1=probY[,1], M2=probY[,2])
rocPlot <- ggplot(df.plot, aes(d=D, m=M2)) +
    geom_roc(n.cuts=50, labels=FALSE) +
    style_roc(theme=theme_grey)
rocPlot <- rocPlot + annotate("text", x=0.75, y=0.25, label=paste("AUC =", round((calc_auc(rocPlot))$AUC, 4)))
print(rocPlot)

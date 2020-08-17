#name: miceImpl
#description: uses mice alrorithm to perform multiple imputation
#language: r
#tags: template, demo
#input: dataframe data [Input data table]
#input: int maxIter = 20 [ scalar giving the number of impute iterations]
#output: dataframe imputedDF [imputed dataset]

require(mice)

imputedData <- mice::mice(data, m = 1,maxit = maxIter)
imputedDF <- mice::complete(imputedData,1)


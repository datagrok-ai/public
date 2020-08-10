#name: miceImpl
#description: uses mice alrorithm to perform multiple imputation
#language: r
#tags: template, demo
#input: dataframe data [Input data table]
#input: int impNum = 5 [number of imputation operations]
#input: int maxIter = 20 [ scalar giving the number of iterations]
#input: int whichImp = 1 [which of the m completed datasets to return]
#output: dataframe outputDF [imputed dataset]

require(mice)

if (whichImp > impNum) {
  whichImp <- 1
}
imputed_Data <- mice::mice(data, m = impNum,maxit = maxIter)
outputDF <- mice::complete(imputed_Data,whichImp)


#name: imputePCAImpl
#description: Impute the missing values of a mixed dataset using the Multiple Correspondence Analysis method
#language: r
#tags: template, demo
#input: dataframe data {slice:false} [Input data table]
#input: column_list columns
#input: int ncp = 9 [number of components used to predict the missing entries]
#input: string method = 'Regularized' [reconstruction formulae]
#input: double regCoeff = 1 [regularization coeff]
#output: dataframe imputedDF [imputed dataset]

require(missMDA)
require(gdata)

# convert all variables to numeric
vars_non_num <- names(data)[!sapply(data, is.numeric)]
bigMap <- mapLevels(data[,c(vars_non_num)])
if (length(vars_non_num) != 0) {
  data <- as.data.frame(sapply(data, as.integer)) }

data <- lapply(data , factor)
imputedDF <- imputePCA(data, ncp = ncp, method = method, coeff.ridge = regCoeff)
imputedDF <- as.data.frame(sapply(imputedDF$completeObs, as.numeric))

mapLevels(imputedDF[,c(vars_non_num)]) <- bigMap
imputedDF <- imputedDF[,columns]
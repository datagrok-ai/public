#name: knnImpl
#description: uses mice alrorithm to perform multiple imputation
#language: r
#tags: template, demo
#input: dataframe data [Input data table]
#input: column_list columns
#input: int k = 10 [number of nearest neighbours]
#output: dataframe imputedDF [imputed dataset]

require(VIM)
require(gdata)

# convert all variables to numeric
vars_non_num <- names(data)[!sapply(data, is.numeric)]
bigMap <- mapLevels(data[,c(vars_non_num)])
if (length(vars_non_num) != 0) {
  data <- as.data.frame(sapply(data, as.integer)) }

Xcolnames <- colnames(data)
imputedDF <- VIM::kNN(data = data, variable = Xcolnames, k = k, trace = F, imp_var = F, numFun = median)

mapLevels(imputedDF[,c(vars_non_num)]) <- bigMap
imputedDF <- imputedDF[,columns]


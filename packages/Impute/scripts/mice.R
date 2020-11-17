#name: miceImpl
#description: uses mice alrorithm to perform multiple imputation
#language: r
#tags: template, demo
#input: dataframe data [Input data table]
#input: column_list columns
#input: int maxIter = 20 [ scalar giving the number of impute iterations]
#output: dataframe imputedDF [imputed dataset]

require(mice)
require(gdata)
require(dplyr)

# convert all variables to numeric
vars_non_num <- names(data)[!sapply(data, is.numeric)]
bigMap <- mapLevels(data[,c(vars_non_num)])
if (length(vars_non_num) != 0) {
  data <- data %>% mutate_at(c(vars_non_num), as.integer)
}

data<-data[rowSums(is.na(data)) != ncol(data), ]
imputedData <- mice::mice(data, m = 1,maxit = maxIter)
imputedDF <- mice::complete(imputedData,1)

mapLevels(imputedDF[,c(vars_non_num)]) <- bigMap
imputedDF <- imputedDF[,columns]
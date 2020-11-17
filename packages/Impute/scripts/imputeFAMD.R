#name: imputeFAMDImpl
#description: Impute the missing values of a mixed dataset using the PCA method "factorial analysis for mixed data"
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
require(dplyr)

# convert all variables to numeric
vars_non_num <- names(data)[!sapply(data, is.numeric)]
vars_int <- names(data)[sapply(data, is.integer)]
data<-data[rowSums(is.na(data)) != ncol(data), ]

if (length(vars_non_num) != 0) {
  bigMap <- mapLevels(data[,c(vars_non_num)])
  data <- data %>% mutate_at(c(vars_non_num), as.factor)
}

imputedDF <- imputeFAMD(data, ncp = ncp, method = method, coeff.ridge = regCoeff)
imputedDF <- as.data.frame(sapply(imputedDF$completeObs, as.numeric))

if (length(vars_non_num) != 0) {
  imputedDF <- imputedDF %>% mutate_at(c(vars_non_num), as.integer)
  mapLevels(imputedDF[,c(vars_non_num)]) <- bigMap
}

imputedDF <- imputedDF %>% mutate_at(c(vars_int), as.integer)
imputedDF <- imputedDF[,columns]
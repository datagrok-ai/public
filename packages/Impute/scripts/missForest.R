#name: missForestImpl
#description: Nonparametric Missing Value Imputation using Random Forest
#language: r
#tags: template, demo
#input: dataframe data
#input: column_list columns
#input: int maxiter = 10
#input: int ntree = 100
#input: bool decreasing = TRUE
#input: bool replace = TRUE
#input: string parallelize  = 'no' {choices : ["no", "variables", "forests"]}
#output: dataframe imputedDF [processed dataframe]

require(missForest)
require(gdata)

# convert all variables to numeric
vars_non_num <- names(data)[!sapply(data, is.numeric)]
bigMap <- mapLevels(data[,c(vars_non_num)])
if (length(vars_non_num) != 0) {
  data <- as.data.frame(sapply(data, as.integer)) }

results <- missForest::missForest(xmis = data,
                                  maxiter = 10,
                                  ntree = 100,
                                  decreasing =,
                                  replace = T,
                                  parallelize = 'no')
imputedDF <- results$ximp

mapLevels(imputedDF[,c(vars_non_num)]) <- bigMap
imputedDF <- imputedDF[,columns]
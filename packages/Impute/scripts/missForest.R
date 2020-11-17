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
require(dplyr)

vars_non_num <- names(data)[!sapply(data, is.numeric)]
vars_int <- names(data)[sapply(data, is.integer)]

data<-data[rowSums(is.na(data)) != ncol(data), ]

if (length(vars_non_num) != 0) { data <- data %>% mutate_at(c(vars_non_num), as.factor) }

results <- missForest::missForest(xmis = data,
                                  maxiter = maxiter,
                                  ntree = ntree,
                                  decreasing =decreasing,
                                  replace = replace,
                                  parallelize = parallelize)

imputedDF <- results$ximp
imputedDF <- imputedDF %>% mutate_at(c(vars_int), as.integer)
imputedDF <- imputedDF[,columns]
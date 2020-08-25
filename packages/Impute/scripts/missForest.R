#name: missForestImpl
#description: Nonparametric Missing Value Imputation using Random Forest
#language: r
#tags: template, demo
#input: dataframe data
#input: int maxiter = 10
#input: int ntree = 100
#input: bool decreasing = TRUE
#input: bool replace = TRUE
#input: string parallelize  = 'no' {choices : ["no", "variables", "forests"]}
#output: dataframe imputedDF [processed dataframe]

require(missForest)

results <- missForest::missForest(xmis = data,
                                  maxiter = 10,
                                  ntree = 100,
                                  decreasing =,
                                  replace = T,
                                  parallelize = 'no')
imputedDF <- results$ximp
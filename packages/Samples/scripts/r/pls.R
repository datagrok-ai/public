#name: PLS
#description: Partial Least Squares (PLS)
#reference: https://en.wikipedia.org/wiki/Partial_least_squares_regression
#language: r
#tags: demo, analysis, pls
#sample: cars.csv
#input: dataframe table [Input data table]
#input: column predict {type:numerical;allowNulls:false} [Predict column]
#input: column_list features {type:numerical;allowNulls:false} [Features column list]
#input: int components = 2 [Number of components]
#output: dataframe predicted {action:join(table)} [Predicted column]
#output: graphics corrLoads

require(plsdepot)

pls1 = plsreg1(table[, features], table[, predict], comps=components)

predicted = data.frame(pls1$y.pred)

corrLoads = plot(pls1, xlab='t1', ylab='t2', main='Correlation Loadings')
print(corrLoads)

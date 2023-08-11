#name: MANOVA
#description: MANOVA with 3 Dependent Variables ((x, y) ~ z)
#reference: https://en.wikipedia.org/wiki/Multivariate_analysis_of_variance
#language: r
#tags: demo, analysis, manova
#sample: beer.csv
#input: dataframe table [Input data table]
#input: column variable1 {type:numerical} [First variable column name]
#input: column variable2 {type:numerical} [Second variable column name]
#input: column variable3 {type:numerical} [Third variable column name]
#output: dataframe manovaResult

require(broom)

assign(paste(variable3, sep=""), table[[variable3]])
formula = as.formula(paste("cbind(table[[variable1]], table[[variable2]]) ~ ", variable3))

manovaResult = manova(formula, table)
manovaResult = tidy(manovaResult)

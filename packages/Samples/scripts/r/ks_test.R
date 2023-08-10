#name: KS Test
#description: Kolmogorov–Smirnov test (K–S test or KS test)
#reference: https://en.wikipedia.org/wiki/Kolmogorov–Smirnov_test
#language: r
#tags: demo
#sample: TSLA.csv
#input: dataframe data [Input data table]
#input: column x {type:numerical} [X axis column name]
#input: column y {type:numerical} [Y axis column name]
#output: double pValue [P-value of ks-statistics]

require(stats)

kstest = ks.test(data[[x]], data[[y]])
pValue = kstest$p.value

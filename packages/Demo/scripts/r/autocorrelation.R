#name: ACF
#description: Autocorrelation function (ACF) plots
#reference: https://en.wikipedia.org/wiki/Autocorrelation
#language: r
#tags: demo, viewers
#sample: TSLA.csv
#input: dataframe data {columns:numerical} [Input data table]
#input: column_list columns {type:numerical; allowNulls:false} [Input data table columns]
#output: graphics acf [Autocorrelation plot]
#test: ACF(ApiTests:getDT(), ['height', 'weight']) != null

require(graphics)

plot = acf(data[columns])

print(plot)

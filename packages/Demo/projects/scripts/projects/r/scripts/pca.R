#name: PCA R
#description: Principal Component Analysis
#reference: https://en.wikipedia.org/wiki/Principal_component_analysis
#language: r
#tags: demo
#sample: cars.csv
#input: dataframe T {columns:numerical} [Input data]
#input: column_list columns {type:numerical;allowNulls:false} [Input data columns]
#input: int numComp = 2 [Number of components]
#input: bool center = TRUE [Indicating whether the variables should be shifted to be zero centered]
#input: bool scale = TRUE [Indicating whether the variables should be scaled to have unit variance before the analysis takes place]
#output: dataframe result {action:join(T)} [PCA components]

require(stats)

# Preform PCA and select required number of components
T.pca <- prcomp(T, center=center, scale.=scale)
result <- predict(T.pca, T)
result = result[,1:numComp]

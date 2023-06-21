#name: PCA Julia
#description: Principal Component Analysis
#reference: https://en.wikipedia.org/wiki/Principal_component_analysis
#language: julia
#tags: demo, hide-suggestions
#sample: cars.csv
#input: dataframe table {columns:numerical} [Input data table]
#input: column_list features {type:numerical;allowNulls:false} [Features columns]
#input: int components = 2 [Number of components]
#output: dataframe result {action:join(table)} [PCA components]

using MultivariateStats

table = Float64.(convert(Array, table))'
model = fit(PCA, table; maxoutdim=components, pratio=0.99999)
result = transform(model, table)
result = convert(DataFrame, result')

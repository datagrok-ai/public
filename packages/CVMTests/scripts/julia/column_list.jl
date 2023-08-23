#name: Julia Column List
#description: column list input
#language: julia
#input: dataframe df
#input: column_list cols
#output: dataframe result

result = select!(df, Not(first(cols)))

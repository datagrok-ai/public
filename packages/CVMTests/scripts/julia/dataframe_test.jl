#name: Julia Dataframe
#description: df input/output
#language: julia
#input: dataframe df
#input: dataframe dfNumerical {columns: numerical}
#input: dataframe dfCategorical {categorical: categorical}
#output: dataframe resultDf
#output: dataframe resultNumerical
#output: dataframe resultCategorical

resultDf = df
resultNumerical = dfNumerical
resultCategorical = dfCategorical
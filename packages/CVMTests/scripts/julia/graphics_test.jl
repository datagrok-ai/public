#name: Julia Graphics
#description: graphics output column input
#language: julia
#input: dataframe df
#input: column xName {type:numerical; allowNulls:false} [Column x]
#input: column yName {type:numerical; allowNulls:false} [Column y]
#output: graphics scatter [scatter plot]

using Plots

p = plot(df[!, xName], df[!, yName])
display(p)

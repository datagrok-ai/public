#name: Julia Int Column
#language: julia
#output: dataframe resultInBound
#output: dataframe resultOutBound
resultInBound = DataFrame((col1=[1000, 10000, 100000]))
resultOutBound = DataFrame((col1=[1000, 10000, typemax(Int32) + 1]))
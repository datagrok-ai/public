#name: R Params Test
#language: r
#tags: test
#input: int i
#input: double d 
#input: bool b
#input: string s
#input: datetime dt
#input_: map m
#input: dataframe df
#input_: column col
#output: int ri
#output: double rd
#output: bool rb
#output: string rs
#output: datetime rdt
#output_: map rm
#output: dataframe rdf
#output_: column rcol
ri = i / 2
rd = d + 60
rb = !b
rs = paste(s, s, sep = '')
rdt = as.Date(dt) + 10
rdf = cbind(df, c(df["height"]))

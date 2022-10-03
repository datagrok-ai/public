#name: R Params Test
#language: r
#tags: test
#input: int i
#input: double d
#input: bool b
#input: string s
#input: datetime dt
#input: dataframe df
#input: column col
#output: int ri
#output: double rd
#output: bool rb
#output: string rs
#output: datetime rdt
#output: dataframe rdf
ri = i / 2
rd = d + 60
rb = !b
rs = paste(s, s, sep = '')
rdt = as.POSIXlt(dt) - 10*86400
rdf = df[c(col)]
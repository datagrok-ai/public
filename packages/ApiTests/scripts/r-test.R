#name: R Params Test
#language: r
#tags: test
#input: int i = 10
#input: double d = -20.1
#input: bool b = false
#input: string s = 'abc'
#input_: datetime dt
#input_: dataframe df
#input_: column col
#output: int ri
#output: double rd
#output: bool rb
#output: string rs
#output_: datetime rdt
#output_: dataframe rdf
#test: ApiTests:getOutput('RParamsTest', 'ri').val == 5
#test: ApiTests:getOutput('RParamsTest', 'rd').val == 39.9
#test: ApiTests:getOutput('RParamsTest', 'rb').val == true
#test: ApiTests:getOutput('RParamsTest', 'rs').val == 'abcabc'

ri = i / 2
rd = d + 60
rb = !b
rs = paste(s, s, sep = '')
#rdt = as.POSIXlt(dt) - 10*86400
#rdf = df[c(col)]
#name: Julia Params Test
#language: julia
#tags: test
#input: int i = 10
#input: double d = -20.1
#input: bool b = false
#input: string s = 'abc'
#input_: datetime dt = '1992-09-20 00:00:00'
#input_: dataframe df {optional: true}
#input_: column col {optional: true}
#output: int ri
#output: double rd
#output: bool rb
#output: string rs
#output_: datetime rdt
#output_: dataframe rdf
#test: ApiTests:getOutput('JuliaParamsTest', 'ri').val == 5
#test: ApiTests:getOutput('JuliaParamsTest', 'rd').val == 39.9
#test: ApiTests:getOutput('JuliaParamsTest', 'rb').val == true
#test: ApiTests:getOutput('JuliaParamsTest', 'rs').val == 'abcabc'

ri = i / 2
rd = d + 60
rb = !b
rs = s^2
#rdf = df[col]
#rdt = dt - timedelta(days=10)

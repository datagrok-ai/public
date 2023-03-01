#name: Python Params Test
#language: python
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
#test: ApiTests:getOutput('PythonParamsTest', 'ri').val == 5
#test: ApiTests:getOutput('PythonParamsTest', 'rd').val == 39.9
#test: ApiTests:getOutput('PythonParamsTest', 'rb').val == true
#test: ApiTests:getOutput('PythonParamsTest', 'rs').val == 'abcabc'

ri = i / 2
rd = d + 60
rb = not b
rs = s + s
#rdf = df[col]
#rdt = dt - timedelta(days=10)

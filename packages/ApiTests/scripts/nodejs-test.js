//name: NodeJS Params Test
//language: nodejs
//tags: test
//input: int i = 10
//input: double d = -20.1
//input: bool b = false
//input: string s = 'abc'
//input_: datetime dt = '1992-09-20 00:00:00'
//input_: dataframe df {optional: true}
//input_: column col {optional: true}
//output: int ri
//output: double rd
//output: bool rb
//output: string rs
//output_: datetime rdt
//output_: dataframe rdf
//test: ApiTests:getOutput('NodeJSParamsTest', 'ri').val == 5
//test: ApiTests:getOutput('NodeJSParamsTest', 'rd').val == 39.90
//test: ApiTests:getOutput('NodeJSParamsTest', 'rb').val == true
//test: ApiTests:getOutput('NodeJSParamsTest', 'rs').val == 'abcabc'

ri = i / 2;
rd = d + 60;
rb = !b;
rs = s + s;
//rdf = df[col]
//rdt = dt - timedelta(days=10)

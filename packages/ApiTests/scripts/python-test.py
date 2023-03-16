#name: Python Params Test
#language: python
#tags: test
#input: int i
#input: double d 
#input: bool b
#input: string s
#input: datetime dt
#input: map m
#input: dataframe df
#input_: column col
#output: int ri
#output: double rd
#output: bool rb
#output: string rs
#output: datetime rdt
#output: map rm
#output: dataframe rdf
#output_: column rcol
ri = i / 2
rd = d + 60
rb = not b
rs = s + s
rdt = dt + timedelta(days=10)
rm = m
rm['b'] = ri
rdf = df
rdf['column'] = rdf['height']
#rdf['column'] = col
#rcol = df['height']

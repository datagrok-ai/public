#name: Python Params Test
#language: python
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
#output_: column rcol
ri = i / 2
rd = d + 60
rb = not b
rs = s + s
rdt = dt + timedelta(days=10)
rdf = df[col]

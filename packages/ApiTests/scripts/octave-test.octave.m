#name: Octave Params Test
#language: octave
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
rb = ~b
rs = strcat(s, s)
rdt = datestr(addtodate(dt, 10, 'days'), 'yyyy-mm-dd')
rdf = df[col]

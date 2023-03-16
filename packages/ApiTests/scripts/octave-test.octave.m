#name: Octave Params Test
#language: octave
#tags: test
#input: int i
#input: double d 
#input: bool b
#input: string s
#input_: datetime dt
#input_: map m
#input: dataframe df
#input_: column col
#output: int ri
#output: double rd
#output: bool rb
#output: string rs
#output_: datetime rdt
#output_: map rm
#output: dataframe rdf
#output_: column rcol
ri = i / 2
rd = d + 60
rb = ~b
rs = strcat(s, s)
#rdt = datestr(addtodate(dt, 10, 'days'), 'yyyy-mm-dd')
rdf = [df, df(:, 1)]

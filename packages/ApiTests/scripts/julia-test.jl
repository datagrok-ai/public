#name: Julia Params Test
#language: julia
#tags: test
#input: int i
#input: double d 
#input_: bool b
#input: string s
#input: datetime dt
#input_: map m
#input: dataframe df
#input_: column col
#output: int ri
#output: double rd
#output_: bool rb
#output: string rs
#output: datetime rdt
#output_: map rm
#output: dataframe rdf
#output_: column rcol
ri = i / 2
rd = d + 60
#rb = !b
rs = s^2
rdt = DateTime(dt, "yyyy-mm-ddTH:M:S.s") + Dates.Day(10)
df.column = col
rdf = df
#rcol = df[:, "height"]
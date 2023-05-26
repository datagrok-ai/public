#name: Julia Params Test
#language: julia
#tags: test
#input: int i
#input_: double d 
#input_: bool b
#input_: string s
#input_: datetime dt
#input_: map m
#input_: dataframe df
#input_: column col
#output: int ri
#output_: double rd
#output_: bool rb
#output_: string rs
#output_: datetime rdt
#output_: map rm
#output_: dataframe rdf
#output_: column rcol
ri = trunc(Int, i / 2)
#rd = d + 60
#rb = !b
#rs = s^2
#rdt = DateTime(dt, "yyyy-mm-ddTH:M:S.s") + Dates.Day(10)
#df.column = col
#rdf = df
#rcol = df[:, "height"]

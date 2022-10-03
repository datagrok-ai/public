#name: Octave Params Test
#language: octave
#tags: test
#input: int i
#input: double d
#input: bool b
#input: string s 
#input: dataframe df [Data table]
#input: column col
#input: column_list col_list
#output: int ri
#output: double rd
#output: bool rb
#output: string rs
#output: dataframe rdf
rb = ~b
rs = [s "-" col_list{1}]
rd = d + 30
ri = rows(df) * columns(df) + i
rdf = {"first_col", "second_col"; 1, 'two'}
#name: Octave Params Test
#language: octave
#tags: test
#input: int i = 10
#input: double d = -20.1
#input: bool b = false
#input: string s = 'abc'
#input_: dataframe df [Data table]
#input_: column col
#input_: column_list col_list
#output: int ri
#output: double rd
#output: bool rb
#output: string rs
#output_: dataframe rdf
#test: ApiTests:getOutput('OctaveParamsTest', 'ri').val == 5
#test: ApiTests:getOutput('OctaveParamsTest', 'rd').val == 39.9
#test: ApiTests:getOutput('OctaveParamsTest', 'rb').val == true
#test: ApiTests:getOutput('OctaveParamsTest', 'rs').val == 'abcabc'

ri = i / 2
rd = d + 60
rb = ~b
rs = strcat(s, s)
#rdf = {"first_col", "second_col"; 1, 'two'}
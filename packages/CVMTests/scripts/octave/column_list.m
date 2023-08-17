#name: Octave Column List
#description: column list input
#language: octave
#input: dataframe df
#input: column_list cols
#output: dataframe result

S = cell2struct(df(1:end,:).', df(1,:) );
S = rmfield(S,cols{1})
result = struct2cell(S)

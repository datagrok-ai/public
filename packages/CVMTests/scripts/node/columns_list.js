//name: NodeJS Column List
//description: column list input
//language: nodejs
//input: dataframe df
//input: column_list cols
//output: dataframe result

result = df.select(...cols.slice(1));

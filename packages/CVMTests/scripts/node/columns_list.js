//name: NodeJS Column List
//description: column list input
//language: nodejs
//input: dataframe df
//input: column_list cols
//output: dataframe result

result = DG.DataFrame.fromColumns(cols.slice(1).map((name) => df.col(name)));

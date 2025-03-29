//name: JavaScript Column List
//description: column list input
//language: javascript
//input: dataframe df
//input: column_list cols
//output: dataframe result

columns = df.columns.remove(cols.toList()[0]);
result = DG.DataFrame.fromColumns(columns.toList());

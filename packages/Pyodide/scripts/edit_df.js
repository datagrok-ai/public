//name: EditDFJS
//description: dataframe input/output
//language: javascript
//input: dataframe input_df
//output: dataframe output_df

const res = input_df.clone();
res.columns.addNewFloat('new_col');
for (let idx = 0; idx < res.rowCount; idx++) {
    res.set('new_col', idx, 1.0);
}
output_df = res;

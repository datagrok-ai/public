#name: EditDFPy
#description: dataframe input/output
#language: pyodide
#input: dataframe input_df
#output: dataframe output_df

rowCount = len(input_df.index)
lst = [1] * rowCount
input_df['new_col'] = lst
output_df = input_df

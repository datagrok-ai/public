#name: Pyodide Column List
#description: column list input
#language: pyodide
#sample: cars.csv
#input: dataframe df
#input: column_list cols
#output: dataframe result
#test: ApiTests:expectTable(PyodideColumnList(OpenFile('System:AppData/Pyodide/cars.d42'), ['model', 'diesel', 'turbo']), OpenFile('System:AppData/Pyodide/cars_2_cols.d42')) //cat: Types

result = df[cols[1:]]

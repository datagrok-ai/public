#name: Pyodide Single Df
#description: df performance
#language: pyodide
#input: dataframe df
#output: dataframe result
#test: ApiTests:expectTable(PyodideSingleDf(OpenFile('System:AppData/Pyodide/cars.d42')), OpenFile('System:AppData/Pyodide/cars.d42')) //cat: Types
result = df
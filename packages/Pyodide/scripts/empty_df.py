#name: Pyodide Empty Dataframe
#language: pyodide
#output: dataframe resultDf
#test: expectTable(PyodideEmptyDataframe(), OpenFile('System:AppData/Pyodide/empty.d42')) //cat: Types
resultDf = pd.DataFrame()

//language: javascript
//name: Pyodide Calc Column Call
//output: dataframe res
//test: ApiTests:expectTable(PyodideCalcColumnCall(), OpenFile('System:AppData/Pyodide/calc_res.d42')) //cat: Types
res = DG.DataFrame.fromCsv('x,y\n1,2\n3,4\n5,6');
await res.columns.addNewCalculated('new', `Pyodide:PyodideCalcColumn(\${x} + \${y})`);

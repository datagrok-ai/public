#name: DataFrameInputCallPy
#language: pyodide
#input: dataframe input_df
#output: dataframe output_df

from dg import execFuncCall
from pyodide.ffi import run_sync

res = run_sync(execFuncCall("Pyodide:EditDFPy", {"input_df": input_df}))

output_df = res["output_df"]

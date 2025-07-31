#name: DatetimeInputCallPy
#language: pyodide
#input: datetime input_datetime
#output: datetime output_datetime

from dg import execFuncCall
from pyodide.ffi import run_sync

res = run_sync(execFuncCall("Pyodide:DatetimeTestPy", {"input_datetime": input_datetime}))

output_datetime = res["output_datetime"]

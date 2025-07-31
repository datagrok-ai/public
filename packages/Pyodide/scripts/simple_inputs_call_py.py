#name: SimpleInputsCallPy
#language: pyodide
#input: bool in1
#input: int in2
#input: double in3
#input: string in4
#output: bool out1
#output: int out2
#output: double out3
#output: string out4

from dg import execFuncCall
from pyodide.ffi import run_sync

res = run_sync(execFuncCall("Pyodide:SimpleInputsPy", {"in1": in1, "in2": in2, "in3": in3, "in4": in4}))
out1 = res["out1"]
out2 = res["out2"]
out3 = res["out3"]
out4 = res["out4"]

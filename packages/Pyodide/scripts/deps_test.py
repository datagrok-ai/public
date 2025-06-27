#name: Pyodide Deps Test
#language: pyodide
#input: bool bool_input = true
#output: bool bool_output
#test: expect(PyodideDepsTest(true), false) //cat: Dependencies
#meta.dependencies: ["openpyxl"]

import openpyxl

bool_output = not bool_input

#name: Pyodide Date
#description: datetime input/output
#language: pyodide
#input: datetime input_datetime
#output: datetime output_datetime
#test: ApiTests:expectDate(PyodideDate(DateTime(1996, 8, 26, 10, 0, 0, 0)), DateTime(1996, 8, 27, 10, 0, 0, 0)) //cat: Types
from datetime import timedelta

output_datetime = input_datetime + timedelta(days=1)


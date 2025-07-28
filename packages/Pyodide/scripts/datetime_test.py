#name: DatetimeTestPy
#description: datetime input/output
#language: pyodide
#input: datetime input_datetime
#output: datetime output_datetime
from datetime import timedelta

output_datetime = input_datetime + timedelta(days=1)

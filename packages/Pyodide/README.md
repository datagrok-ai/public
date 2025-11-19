# Pyodide

The Pyodide package integrates [Pyodide](https://pyodide.org/), enabling the execution of Python scripts directly in the browser.

## Features

- **In-browser Python execution**: Run Python code directly within Datagrok without the need for external servers.
- **Seamless integration**: Easily integrate Python functions into Datagrok workflows and visualizations.
- **Data interoperability**: Pass data between Datagrok and Python functions effortlessly.
- **Asynchronous execution**: Execute Python functions asynchronously to maintain responsive UIs.

All you need to do is to set the language to "pyodide", like in this example:

```python
#name: PyodideTest
#description: dataframe input/output
#language: pyodide
#input: dataframe input_df
#output: dataframe output_df

rowCount = len(input_df.index)
lst = [1] * rowCount
input_df['new_col'] = lst
output_df = input_df
```

See also:
* [Scripting](https://datagrok.ai/help/compute/scripting/)

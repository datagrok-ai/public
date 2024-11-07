#name: Pyodide Bool
#language: pyodide
#input: bool bool_input = true
#output: bool bool_output
#test: expect(PyodideBool(true), false) //cat: Types

bool_output = not bool_input

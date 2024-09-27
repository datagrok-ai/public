#name: Pyodide String
#language: pyodide
#input: string string_input
#output: string string_output
#test: expect(PyodideString('Hello world'), 'Hello world!') //cat: Types

string_output = string_input + '!'

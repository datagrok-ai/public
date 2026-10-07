#name: Inline Env

#language: pyodide
#environment: channels: [conda-forge], dependencies: [python=3.12, {pip: [cobra]}]
#input: int x
#output: int y
y = x

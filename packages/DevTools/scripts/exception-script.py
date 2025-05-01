#name: ExceptionScriptPython
#language: python
#input: int a
#output: int out

if a == 0:
    raise Exception('exception')
else:
    a += 1
out = a

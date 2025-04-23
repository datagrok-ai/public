//name: ExceptionScriptNode
//language: nodejs
//input: int a
//output: int out

if (a == 0) // use type coercion
  throw 'exception';
else
  a++;
out = a;

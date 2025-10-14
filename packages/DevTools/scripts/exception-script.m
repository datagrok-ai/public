#name: ExceptionScriptOctave
#language: octave
#input: int a
#output: int out

if (a == 0)
  error("exception");
else
  a += 1;
endif
out = a;

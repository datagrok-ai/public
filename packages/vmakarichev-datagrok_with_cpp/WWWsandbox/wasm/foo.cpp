#include <emscripten.h>

extern "C" {
  int foo(int n);
}

//name: fib
//input: int n
//output: double res [_callResult]
EMSCRIPTEN_KEEPALIVE
int foo(int n)
{    
    return n * 2;
}
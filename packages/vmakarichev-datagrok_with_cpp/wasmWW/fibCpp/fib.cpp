#include <emscripten.h>

extern "C" {
  int fib(int n);
}

//name: fib
//input: int n
//output: double res [_callResult]
EMSCRIPTEN_KEEPALIVE
int fib(int n)
{
    if(n <= 1)
      return 1;
    return fib(n - 1) + fib(n - 2);
}
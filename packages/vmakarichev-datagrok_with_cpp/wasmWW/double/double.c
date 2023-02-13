#include <emscripten.h>

EMSCRIPTEN_KEEPALIVE
void doubleArray(int * arr, int len)
{
    for(int i = 0; i < len; i++)
      arr[i] *= 2;
}
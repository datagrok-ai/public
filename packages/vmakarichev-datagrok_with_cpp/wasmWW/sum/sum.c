#include <emscripten.h>

EMSCRIPTEN_KEEPALIVE
int sum(int * arr, int len)
{
    int res = 0;

    for(int i = 0; i < len; i++)
      res += arr[i];
    
    return res;
}
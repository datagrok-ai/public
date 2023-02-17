#include <emscripten.h>

extern "C" {
  int fillArray(float * arr, int arrLength);
}

EMSCRIPTEN_KEEPALIVE
int fillArray(float * arr, int arrLength)
{
    for(int i = 0; i < arrLength; i++)
      arr[i] = 1.0f / (i + 1);
      
    return 0;
}
#include <emscripten.h>
extern "C" {
    float minOfArr(float * arr, int arrLength);
    float maxOfArr(float * arr, int arrLength);
}


#include "stat/statistics.h"

//name: minOfColumn
//input: dataframe df
//input: column col
//output: double num
EMSCRIPTEN_KEEPALIVE
float minOfArr(float * arr, int arrLength)
{
    return sta::minOfArray(arr, arrLength);
}

//name: maxOfColumn
//input: dataframe df
//input: column col
//output: double num
EMSCRIPTEN_KEEPALIVE
float maxOfArr(float * arr, int arrLength)
{
    return sta::maxOfArray(arr, arrLength);
}
#include <emscripten.h>
extern "C" {
    float minOfArrOfFloats(float * arr, int arrLength);
    float maxOfArrOfFloats(float * arr, int arrLength);
    int minOfArrOfInts(int * arr, int arrLength);
    int maxOfArrOfInts(int * arr, int arrLength);
}


#include "stat/indicators.h"

//name: minOfColumnOfDoubles
//input: dataframe df
//input: column col
//output: double num
EMSCRIPTEN_KEEPALIVE
float minOfArrOfFloats(float * arr, int arrLength)
{
    return ind::minOfArray(arr, arrLength);
}

//name: maxOfColumnOfDoubles
//input: dataframe df
//input: column col
//output: double num
EMSCRIPTEN_KEEPALIVE
float maxOfArrOfFloats(float * arr, int arrLength)
{
    return ind::maxOfArray(arr, arrLength);
}

//name: minOfColumnOfInts
//input: dataframe df
//input: column col
//output: int num
EMSCRIPTEN_KEEPALIVE
int minOfArrOfInts(int * arr, int arrLength)
{
    return ind::minOfArray(arr, arrLength);
}

//name: maxOfColumnOfInts
//input: dataframe df
//input: column col
//output: int num
EMSCRIPTEN_KEEPALIVE
int maxOfArrOfInts(int * arr, int arrLength)
{
    return ind::maxOfArray(arr, arrLength);
}
#include <emscripten.h>

EMSCRIPTEN_KEEPALIVE
void doubleArray(float * array, int length)
{
    for(int i = 0; i < length; i++)
        array[i] = array[i] * 2;

} // doubleArray

EMSCRIPTEN_KEEPALIVE
int minOfArray(int * array, int length)
{
    int result = array[0];

    for(int i = 1; i < length; i++)
        if( result > array[i] )
            result = array[i];
    
    return result;
} // minOfArray


EMSCRIPTEN_KEEPALIVE
void sumOfArrays(float * arr1, float * arr2, float * sum, int length)
{
    for(int i = 0; i < length; i++)
        sum[i] = arr1[i] + arr2[i];
} // sumOfArrays
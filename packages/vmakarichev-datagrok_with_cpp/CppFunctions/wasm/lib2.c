#include <emscripten.h>

//name: doubleColumn
//input: column col
//output: column result
//map: array is col
//map: length is col.length
//update: result is array
EMSCRIPTEN_KEEPALIVE
void doubleArray(float * array, int length)
{
    for(int i = 0; i < length; i++)
        array[i] = array[i] * 2;

} // doubleArray

//name: sumOfColumns
//input: column col1
//input: column col2
//output: column result
//map: arr1 is col1
//map: arr2 is col2
//map: sum is new(length)
//map: length is col1.length
//update: result is sum
EMSCRIPTEN_KEEPALIVE
void sumOfArrays(float * arr1, float * arr2, float * sum, int length)
{
    for(int i = 0; i < length; i++)
        sum[i] = arr1[i] + arr2[i];
} // sumOfArrays
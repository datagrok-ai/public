#include <emscripten.h>

//name: sumOfIntegers
//input: int x
//input: int y
//output: int z
//map: a is x
//map: b is y
EMSCRIPTEN_KEEPALIVE
int sumOfInts(int a, int b)
{
    return a + b;
} // sumOfInts

//name: minOfColumn
//input: dataframe df
//input: column col
//output: int num
//map: array is col
//map: length is col.length
EMSCRIPTEN_KEEPALIVE
int minOfArray(int * array, int length)
{
    int result = array[0];

    for(int i = 1; i < length; i++)
        if( result > array[i] )
            result = array[i];
    
    return result;
} // minOfArray


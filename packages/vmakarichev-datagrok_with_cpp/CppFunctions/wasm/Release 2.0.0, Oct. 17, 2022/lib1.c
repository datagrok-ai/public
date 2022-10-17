#include <emscripten.h>

//name: sumOfIntegers
//input: int x
//input: int y
//output: int z
EMSCRIPTEN_KEEPALIVE
int sumOfInts(int a, int b)
{
    return a + b;
} // sumOfInts

//name: minOfColumn
//input: dataframe df
//input: column col
//output: int num
EMSCRIPTEN_KEEPALIVE
int minOfArray(int * array, int arrayLength)
{
    int result = array[0];

    for(int i = 1; i < arrayLength; i++)
        if( result > array[i] )
            result = array[i];
    
    return result;
} // minOfArray


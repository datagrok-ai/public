#include <emscripten.h>

//name: doubleColumn (non-export)
//input: column col
//output: column result [array]
EMSCRIPTEN_KEEPALIVE
void doubleArray(float * array, int arrayLength)
{
    for(int i = 0; i < arrayLength; i++)
        array[i] = array[i] * 2;

}

//name: sumOfColumns (non-export)
//input: column col1
//input: column col2
//output: column result [sum new(col1.length)]
EMSCRIPTEN_KEEPALIVE
void sumOfArrays(float * arr1, int arr1Length, float * arr2, int arr2Length, float * sum, int sumLength)
{
    for(int i = 0; i < sumLength; i++) {
        sum[i] = arr1[i] + arr2[i];
    }
} 

//name: sumAndProdOfColumns (non-export)
//input: column col1
//input: column col2
//output: columns result [sum new(col1.length), prod new(...)]
EMSCRIPTEN_KEEPALIVE
void sumAndProdOfArrays(float * arr1, int arr1Length, float * arr2, int arr2Length, float * sum, int sumLength, float * prod, int prodLength)
{
    for(int i = 0; i < sumLength; i++) {
        sum[i] = arr1[i] + arr2[i];
        prod[i] = arr1[i] * arr2[i];
    }
}

//name: powersOfInts (non-export)
//input: column col
//output: columns result [squared new(col.length), cubic new(...), fourth new(...)]
EMSCRIPTEN_KEEPALIVE
void powers(int * arr, int arrLength, int * squared, int squaredLength, int * cubic, int cubicLength, int * fourth, int fourthLength) {
    for(int i = 0; i < arrLength; i++) {
        squared[i] = arr[i] * arr[i];
        cubic[i] = squared[i] * arr[i];
        fourth[i] = cubic[i] * arr[i];
    }
}
// This file contains C++-functions that are exported to wasm.

// The tool Emscripten is applied (the header emscripten.h is included
// and each exported function is marked by EMSCRIPTEN_KEEPALIVE).

// Also, each function has a special DATAGROK annotation for C++-functions.
// This approach provides further usage of C++-to-wasm export script that 
// performes all routine steps. 

#include <emscripten.h>

#include<cmath>
using namespace std;

// The following provides convenient naming of the exported functions.
extern "C" {
    int sum(int a, int b);

    float maxFloatCol(float * col, int colLength);

    int maxIntCol(int * col, int colLength);

    int addFloatCols(float * col1, int col1Length, 
                     float * col2, int col2Length,
                     float * sum, int sumLength);

    int addIntCols(int * col1, int col1Length, 
                   int * col2, int col2Length,
                   int * sum, int sumLength);

    int doubledInts(int * cols, int colsRowsCount, int colsColsCount,
                    int * res, int resRowsCount, int resColsCount);
                    
    int doubledFloats(float * cols, int colsRowsCount, int colsColsCount,
                      float * res, int resRowsCount, int resColsCount);
};

//name: sum
//input: int a
//input: int b
//output: int res [_callResult]
EMSCRIPTEN_KEEPALIVE  
int sum(int a, int b)
{
    return a + b;
}

//name: maxFloatCol
//input: dataframe df
//input: column col
//output: double max [_callResult]
EMSCRIPTEN_KEEPALIVE
float maxFloatCol(float * col, int colLength)
{
    float res = col[0];

    for(int i = 1; i < colLength; i++)
        res = fmax(res, col[i]);
    
    return res;
}

//name: maxIntCol
//input: dataframe df
//input: column col
//output: int max [_callResult] 
EMSCRIPTEN_KEEPALIVE
int maxIntCol(int * col, int colLength)
{
    int res = col[0];

    for(int i = 1; i < colLength; i++)
        res = fmax(res, col[i]);
    
    return res;
}

//name: addFloatCols
//input: dataframe df
//input: column col1
//input: column col2
//output: column sum [new(col1.rowCount)]
EMSCRIPTEN_KEEPALIVE
int addFloatCols(float * col1, int col1Length, 
                 float * col2, int col2Length,
                 float * sum, int sumLength)
{
    for(int i = 0; i < sumLength; i++)
        sum[i] = col1[i] + col2[i];

    return 0;
}

//name: addIntCols
//input: dataframe df
//input: column col1
//input: column col2
//output: column sum [new(col1.rowCount)]
EMSCRIPTEN_KEEPALIVE
int addIntCols(int * col1, int col1Length, 
               int * col2, int col2Length,
               int * sum, int sumLength)
{
    for(int i = 0; i < sumLength; i++)
        sum[i] = col1[i] + col2[i];

    return 0;
}

//name: doubledInts
//input: dataframe table
//input: column_list cols
//output: column_list components [new(cols.rowCount, cols.columnCount)]
//output: dataframe result [components]
EMSCRIPTEN_KEEPALIVE
int doubledInts(int * cols, int colsRowCount, int colsColCount,
                int * res, int resRowCount, int resColCount)
{
    for(int i = 0; i < colsRowCount * colsColCount; i++)
        res[i] = 2 * cols[i];
    
    return 0;
}

//name: doubledFloats
//input: dataframe table
//input: column_list cols
//output: column_list components [new(cols.rowCount, cols.columnCount)]
//output: dataframe result [components]
EMSCRIPTEN_KEEPALIVE
int doubledFloats(float * cols, int colsRowCount, int colsColCount,
                  float * res, int resRowCount, int resColCount)
{
    for(int i = 0; i < colsRowCount * colsColCount; i++)
        res[i] = 2.0f * cols[i];
    
    return 0;
}
#include <emscripten.h>

extern "C" {
    int basicExample(float t0, float t1, int timesCount, float parameter, int resultCount,
                 float * time, int timeLength,
                 float * result, int resultRowCount, int resultColumnCount);
}

#include "wrapper.h"


//name: example
//input: double t0
//input: double t1
//input: int timesCount
//input: double parameter
//input: int resultCount
//output: column time [new(timesCount)]
//output: column_list result [new(timesCount, resultCount)]
EMSCRIPTEN_KEEPALIVE
int basicExample(float t0, float t1, int timesCount, float parameter, int resultCount,
                 float * time, int timeLength,
                 float * result, int resultRowCount, int resultColumnCount)
{
    if (timesCount < 2)
        return -1;

    double step = (t1 - t0) / (timesCount - 1);

    time[0] = t0;

    for(int i = 1; i < timesCount; i++)    
        time[i] = time[i - 1] + step;

    float params[1];
    params[0] = parameter;

    return solverWrapper(time, timesCount, params, 1, result, resultRowCount, resultColumnCount);
}
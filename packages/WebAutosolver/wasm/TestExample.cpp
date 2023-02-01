// TestExample.cpp

#include <cmath>
using namespace std;

#include <emscripten.h>

extern "C" {
     int solveTestExample(float initial, float final, float step,
        float x, float y, 
        float param1, float param2, 
        float * result, int resultRowCount, int resultColCount) noexcept;
}

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "odes.h"
using namespace odes;
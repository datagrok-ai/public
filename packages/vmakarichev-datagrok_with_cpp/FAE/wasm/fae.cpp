// fae.cpp

#include <emscripten.h>

extern "C" {
    int solveFAE(float t0, float t1, float h, int timesCount, int varsCount,
                 float * result, int resultRowCount, int resultColCount) noexcept;
}

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "odeSolver.h"
using namespace ode;

namespace fae
{
	// overall tollerance
    float tol = 0.00005f; 

	float yInitial[] = { // initial conditions are removed!
	};	
    
	VectorXd f(double _time, VectorXd & y) noexcept
	{
		// original formulas are removed!
	} // f

}; //fae

//name: solveFAE
//input: double t0
//input: double t1
//input: double h
//input: int timesCount
//input: int varsCount
//output: column_list result [new(timesCount, varsCount)]
//output: dataframe solution [result]
EMSCRIPTEN_KEEPALIVE
int solveFAE(float t0, float t1, float h, int timesCount, int varsCount,
  float * result, int resultRowCount, int resultColCount) noexcept
{
    using namespace fae;

    return solveODE(f, t0, t1, h, yInitial, tol, result, timesCount, varsCount);
} // solveFAE


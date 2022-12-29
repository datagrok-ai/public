// stiffExample.cpp

#include <emscripten.h>

extern "C" {
    int stiffExample(float t0, float t1, int timesCount, int varsCount,
	                 float * result, int resultRowCount, int resultColCount);

	int stiffExampleRK(float t0, float t1, int timesCount, int varsCount,
	                   float * result, int resultRowCount, int resultColCount);
}

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "odeSolver.h"
using namespace ode;

namespace stiffex 
{
	VectorXd F(double t, VectorXd& y)
	{
		VectorXd res(y.size());

		res(0) = -5.0 * y(0) + 3.0 * y(1);
		res(1) = 100.0 * y(0) - 301.0 * y(1);

		return res;
	}

	MatrixXd J(double t, VectorXd& y, double eps)
	{
		MatrixXd res(y.size(), y.size());

		res << -5.0, 3.0,
			100.0, -301.0;

		return res;
	}

	VectorXd T(double t, VectorXd& y, double eps)
	{
		VectorXd res(y.size());

		res << 0.0,
			0.0;

		return y;
	}

	VectorXd exactSolution(double t, unsigned dim)
	{
		VectorXd res(dim);

		res(0) = 52.96 * exp(-3.9899 * t) - 0.67 * exp(-302.0101 * t);
		res(1) = 17.83 * exp(-3.9899 * t) + 65.99 * exp(-302.0101 * t);

		return res;
	}
}; // stiffex

//name: stiffExample
//input: double t0
//input: double t1
//input: int timesCount
//input: int varsCount
//output: column_list result [new(timesCount, varsCount)]
//output: dataframe solution [result]
EMSCRIPTEN_KEEPALIVE
int stiffExample(float t0, float t1, int timesCount, int varsCount,
  float * result, int resultRowCount, int resultColCount)
{
	using namespace stiffex;

	if (timesCount < 2)
        return -1;

	resultRowCount = timesCount;
	resultColCount = varsCount;

	float step = (t1 - t0) / (timesCount - 1);

	double* times = new double[resultRowCount];
	double* solution = new double[resultRowCount * (resultColCount - 1)];	

	times[0] = result[0] = t0;
	for (int i = 1; i < resultRowCount; i++)
		times[i] = result[i] = result[i - 1] + step;

	VectorXd y = exactSolution(t0, resultColCount - 1);

	int resultCode = oneStepSolver(F, T, J, times, resultRowCount, y.data(), resultColCount - 1, solution);

	if (resultCode != NO_ERRORS)
		return resultCode;

	for (int j = 1; j < resultColCount; j++)
		for (int i = 0; i < resultRowCount; i++)
			result[i + j * resultRowCount] = solution[i + resultRowCount * (j - 1)];

	delete[] times;
	delete[] solution;

	return resultCode;
} // stiffExample


//name: stiffExampleRK
//input: double t0
//input: double t1
//input: int timesCount
//input: int varsCount
//output: column_list result [new(timesCount, varsCount)]
//output: dataframe solution [result]
EMSCRIPTEN_KEEPALIVE
int stiffExampleRK(float t0, float t1, int timesCount, int varsCount,
  float * result, int resultRowCount, int resultColCount)
{
	using namespace stiffex;

	if (timesCount < 2)
        return -1;

	resultRowCount = timesCount;
	resultColCount = varsCount;

	float step = (t1 - t0) / (timesCount - 1);

	double* times = new double[resultRowCount];
	double* solution = new double[resultRowCount * (resultColCount - 1)];

	times[0] = result[0] = t0;
	for (int i = 1; i < resultRowCount; i++)
		times[i] = result[i] = result[i - 1] + step;

	VectorXd y = exactSolution(t0, resultColCount - 1);

	int resultCode = oneStepSolver(F, times, resultRowCount, y.data(), resultColCount - 1, getNextPointRK4, solution);
		
	if (resultCode != NO_ERRORS)
		return resultCode;

	for (int j = 1; j < resultColCount; j++)
		for (int i = 0; i < resultRowCount; i++)
			result[i + j * resultRowCount] = solution[i + resultRowCount * (j - 1)];

	delete[] times;
	delete[] solution;

	return 0;
} // stiffExampleRK
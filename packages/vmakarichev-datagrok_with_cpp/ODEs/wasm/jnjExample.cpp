// jnjExample.cpp

#include <emscripten.h>

extern "C" {
    int jnjStiff(float t0, float t1, int timesCount, int varsCount,
	             float * result, int resultRowCount, int resultColCount);

	int jnjRK4(float t0, float t1, int timesCount, int varsCount,
	           float * result, int resultRowCount, int resultColCount);
}

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "odeSolver.h"
using namespace ode;

namespace demoJNJ
{
		VectorXd F(double t, VectorXd & y)
		{
			VectorXd res(y.size());

		}

		MatrixXd J(double t, VectorXd & y, double eps)
		{
			MatrixXd res(y.size(), y.size());

			VectorXd val = F(t, y);

			VectorXd yDer = y;

			for (int i = 0; i < y.size(); i++)
			{
				yDer(i) += eps;

				res.col(i) = (F(t, yDer) - val) / eps;

				yDer(i) -= eps;
			}

			return res;
		}

		VectorXd T(double t, VectorXd & y, double eps)
		{
			return (F(t + eps, y) - F(t, y)) / eps;
		}

		VectorXd exactSolution(double t, unsigned dim)
		{
			VectorXd res(dim);

			

			return res;
		}
}; //demoJNJ

//name: jnjStiff
//input: double t0
//input: double t1
//input: int timesCount
//input: int varsCount
//output: column_list result [new(timesCount, varsCount)]
//output: dataframe solution [result]
EMSCRIPTEN_KEEPALIVE
int jnjStiff(float t0, float t1, int timesCount, int varsCount,
  float * result, int resultRowCount, int resultColCount)
{
	using namespace demoJNJ;

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
} // jnjStiff


//name: jnjRK4
//input: double t0
//input: double t1
//input: int timesCount
//input: int varsCount
//output: column_list result [new(timesCount, varsCount)]
//output: dataframe solution [result]
EMSCRIPTEN_KEEPALIVE
int jnjRK4(float t0, float t1, int timesCount, int varsCount,
  float * result, int resultRowCount, int resultColCount)
{
	using namespace demoJNJ;

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
} // jnjRK4
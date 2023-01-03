// stiffExample.cpp

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "odeSolver.h"
using namespace ode;

#include "demos.h"

namespace demoJNJ
{
		VectorXd F(double t, VectorXd & y)
		{
			VectorXd res(y.size());

			// JNJ formulas are removed!
			return res;
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

			// JNJ formulas are removed!

			return res;
		}
}; //demoJNJ

// demo of solving stiff example
int jnjStiff(float t0, float t1, float step,
	float* result, int resultRowCount, int resultColCount)
{
	using namespace demoJNJ;

	double* times = new double[resultRowCount];
	double* solution = new double[resultRowCount * (resultColCount - 1)];

	times[0] = result[0] = t0;
	for (int i = 1; i < resultRowCount; i++)
		times[i] = result[i] = result[i - 1] + step;

	VectorXd y = exactSolution(t0, resultColCount - 1);

	int resultCode = fixedStepSolver(F, T, J, times, resultRowCount, y.data(), resultColCount - 1, solution);

	if (resultCode != NO_ERRORS)
		return resultCode;

	for (int j = 1; j < resultColCount; j++)
		for (int i = 0; i < resultRowCount; i++)
			result[i + j * resultRowCount] = solution[i + resultRowCount * (j - 1)];

	delete[] times;
	delete[] solution;

	return 0;
} // jnjStiff


// demo of solving stiff example
int jnjRK4(float t0, float t1, float step,
	float* result, int resultRowCount, int resultColCount)
{
	using namespace demoJNJ;

	double* times = new double[resultRowCount];
	double* solution = new double[resultRowCount * (resultColCount - 1)];

	times[0] = result[0] = t0;
	for (int i = 1; i < resultRowCount; i++)
		times[i] = result[i] = result[i - 1] + step;

	VectorXd y = exactSolution(t0, resultColCount - 1);

	int resultCode = fixedStepSolver(F, times, resultRowCount, y.data(), resultColCount - 1, getNextPointRK4, solution);

	if (resultCode != NO_ERRORS)
		return resultCode;

	for (int j = 1; j < resultColCount; j++)
		for (int i = 0; i < resultRowCount; i++)
			result[i + j * resultRowCount] = solution[i + resultRowCount * (j - 1)];

	delete[] times;
	delete[] solution;

	return 0;
} // jnjRK4
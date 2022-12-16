// testNewtonSolverND.cpp

#include<iostream>
#include<vector>
using namespace std;

typedef vector<double> Vec;

#include "test.h"

Vec F(Vec & x)
{
	Vec y(x.size());

	/*y[0] = x[0] * x[0] - x[0];
	y[1] = exp(-x[1]) - x[1];*/

	y[0] = x[0] * x[0] + x[1] * x[1] - 5;
	y[1] = x[0] * x[1] - 2;

	return y;
}

Vec getRoots(Vec(*F)(Vec &), Vec & x0, double precision, double eps)
{
	int N = x0.size();
	Vec xPrev(x0);
	Vec xCurrent(N);

	int iter = 1;

	while (true)
	{
		auto res = F(xPrev);			

		for (int i = 0; i < N; i++)
		{
			Vec xDer(xPrev);

			xDer[i] += eps;

			auto resDer = F(xDer);

			double lambda = (resDer[i] - res[i]) / eps;

			//cout << "  " << lambda;

			xCurrent[i] = xPrev[i] - res[i] / lambda;
		}

		cout << "\n xCurrent:\n";
		for (int i = 0; i < N; i++)
			cout << "  " << xCurrent[i];

		double deviation = fabs(xCurrent[0] - xPrev[0]);

		for (int i = 1; i < N; i++)
			deviation = fmax(deviation, fabs(xCurrent[i] - xPrev[i]));

		cout << "\n  iter no. " << iter << "   mad = " << deviation << endl;

		if (deviation < precision)
			break;

		iter++;

		for (int i = 0; i < N; i++)
			xPrev[i] = xCurrent[i];
	}

	return xCurrent;
}

// test of solver of the equation F(x) = 0 using Newton's method: multi-dimensional case
void testNewtonSolverND()
{
	int N = 2;
	double precision = 0.0001;
	double eps = 0.001;

	Vec x0(N);
	x0[0] = 1000000; x0[1] = -20000;

	auto solution = getRoots(F, x0, precision, eps);

	cout << "\nRoots:\n";
	for (int i = 0; i < N; i++)
		cout << "  " << solution[i];
}
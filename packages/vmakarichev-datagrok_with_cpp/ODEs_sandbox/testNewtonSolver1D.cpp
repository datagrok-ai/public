#include <cmath>
#include <iostream>
using namespace std;

#include "test.h"

double g(double x)
{
	//return x * x;

	//return exp(-x);

	return 1.0 / x;

	return cos(x);

	return sin(x);
}

// Returns root of the equation x = g(x) obtained using Newton method
double getRootOneDim(double(*g)(double), double initialPoint, double precision, double eps)
{
	double xPrev = initialPoint;
	double xCurrent = 0.0;
	double val = 0.0;
	int iter = 1;

	while (true)
	{
		val = g(xPrev);

		xCurrent = xPrev - (val - xPrev) * eps / (g(xPrev + eps) - val - eps);

		cout << "\n  iter no. " << iter << "  x = " << xCurrent << "   deviation = " << fabs(xCurrent - xPrev) << endl;

		if (fabs(xCurrent - xPrev) < precision)
			break;

		xPrev = xCurrent;
		iter++;
	}

	return xCurrent;
}

// test of solver of the equation x = g(x) using Newton's method: one-dimensional case
void testNewtonSolver1D()
{
	double x0 = 100;
	double precision = 0.01;
	double eps = 0.0000001;

	double root = getRootOneDim(g, x0, precision, eps);

	cout << "\nInitial point: x0 = " << x0
		<< "\nroot = " << root
		<< "\ng(root) = " << g(root) << endl;
}
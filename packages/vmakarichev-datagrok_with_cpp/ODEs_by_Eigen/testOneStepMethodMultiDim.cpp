// testOneStepMethodMultiDim.cpp

// Implementation of the function testOneStepMethodMultiDim() - test of one-step method: multi-dimensional case

#include<iostream>
#include<iomanip>
using namespace std;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "odeSolver.h"
using namespace ode;

//VectorXd f(double t, VectorXd & y)
//{
//	VectorXd res(y.size());
//
//	res(0) = 2.0 * t;
//	res(1) = 3.0 * t * t;
//
//	return res;
//}
//
//VectorXd exact(double t, unsigned dim)
//{
//	VectorXd res(dim);
//
//	res(0) = t * t; 
//	res(1) = t * t * t;
//
//	return res;
//}

/*VectorXd f(double t, VectorXd & y)
{
	VectorXd res(y.size());

	res(0) = 2.0 * exp(2.0 * t);
	res(1) = -3.0 * exp(-3.0 * t);

	return res;
}

VectorXd exact(double t, unsigned dim)
{
	VectorXd res(dim);

	res(0) = exp(2.0 * t);
	res(1) = exp(-3.0 * t);

	return res;
}*/

VectorXd f(double t, VectorXd & y)
{
	VectorXd res(y.size());

	res(0) = sin(t);
	res(1) = cos(t);

	return res;
}

VectorXd exact(double t, unsigned dim)
{
	VectorXd res(dim);

	res(0) = -cos(t);
	res(1) = sin(t);

	return res;
}

//VectorXd f(double t, VectorXd & y)
//{
//	VectorXd res(y.size());
//
//	res(0) = 2.0 * y(1);
//	res(1) = 2.0 * y(0);
//
//	return res;
//}
//
//VectorXd exact(double t, unsigned dim)
//{
//	VectorXd res(dim);
//
//	res(0) = exp(2.0 * t) + exp(-2.0 * t);
//	res(1) = exp(2.0 * t) - exp(-2.0 * t);
//
//	return res;
//}

// test of one-step method: multi-dimensional case
void testOneStepMethodMultiDim()
{
	cout << "Test of one step method: multi-dimensional case.\n";
	
	/*VectorXd y(2);

	y << 0.0, 0.0;

	double t = 0.0;
	double h = 0.001;	

	int n = 100;

	cout << "\nSolution:\n";
	cout << "  i = " << 0 << "   t = " << t << "\n   y =\n" << y << endl;

	for (int i = 1; i <= n; i++)
	{
		y = getNextPointRK4(f, t, y, h);
		t += h;

		cout << "\n  i = " << i << "   t = " << t << "\n   y = " << y.transpose() << endl;
	}*/

	double h = 0.01;
	const unsigned N = 100001;
	const unsigned DIM = 2;

	double * times = new double[N];
	double * solution = new double[N * DIM];

	times[0] = 0.0;
	VectorXd yInitial = exact(times[0], DIM);

	for (int i = 1; i < N; i++)
		times[i] = times[i - 1] + h;

	cout << "Solving ODE result code: ";
	cout << oneStepSolver(f, times, N, yInitial.data(), DIM, getNextPointRK5, solution) << endl;

	cout << "\n   Time         y1(approx)         y2(approx)         y1(exact)         y2(exact)      Error(absolute)\n";
	for (int i = 0; i < N; i++)
	{
		VectorXd yExact = exact(times[i], DIM);
		cout << setw(7) << setiosflags(ios::right) << times[i]
			<< setw(19) << solution[i]
			<< setw(18) << solution[i + N]
			<< setw(18) << yExact(0)
			<< setw(18) << yExact(1)
			<< setw(22) << fmax(fabs(solution[i] - yExact(0)), fabs(solution[i + N] - yExact(1)))
			<< endl;
	}

	delete[] times;
	delete[] solution;
}
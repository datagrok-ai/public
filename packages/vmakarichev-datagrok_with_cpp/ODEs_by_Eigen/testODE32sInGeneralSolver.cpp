// testODE32sInGeneralSolver.cpp

#include<iostream>
#include<iomanip>
#include<ctime>
#include<fstream>
using namespace std;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "odeSolver.h"
using namespace ode;

#include "tests.h"

namespace testODE32sGSsimple
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
} // testODE32sGSsimple

// test of ode32s solver in generalized solver: multi-dimensional case
void testODE32sInGeneralSolver()
{
	using namespace testODE32sGSsimple;
	
	double t0 = 0.0;
	double t1 = 10000.0;
	double h = 0.0002;
	unsigned n = static_cast<unsigned>((t1 - t0) / h) + 1;
	unsigned d = 2;

	cout << "Test of ODE32s solver in the generalized form.\n";

	cout << "\nThe initial problem is solved on the segment [t0, t1] with the step h,\nwhere\n"
		<< "   t0 = " << t0
		<< "\n   t1 = " << t1
		<< "\n   h = " << h
		<< "\n\n   number of points = " << n
		<< endl;

	double* times = new double[n];
	double* solution = new double[n * d];

	times[0] = t0;
	for (unsigned i = 1; i < n; i++)
		times[i] = times[i - 1] + h;

	VectorXd y = exactSolution(t0, d);	

	cout << "\nSolving system of equations...\n";

	cout << "\nResult code: "
		<< fixedStepSolver(F, T, J, times, n, y.data(), d, solution)
		//<< oneStepSolver(F, times, n, y.data(), d, getNextPointRK4, solution)
		<< endl;

	double mad = 0.0;

	cout << "\n  time:         y1(approx):          y1(exact):         y2(approx):          y2(exact):            max error:\n";

	for (unsigned i = 0; i < n; i++)
	{
		VectorXd y = exactSolution(times[i], d);

		double error = fmax(fabs(y(0) - solution[i]), fabs(y(1) - solution[i + n]));

		mad = fmax(mad, error);

		/*cout << setiosflags(ios::right)
			<< setw(7) << times[i]
			<< setw(20) << solution[i]
			<< setw(20) << y(0)
			<< setw(20) << solution[i + n]			
			<< setw(20) << y(1)
			<< setw(22) << error
			<< endl;*/
	}

	cout << "\n----------------------------------------------------------------------------\n"
		<< "\nMAX ERROR: " << mad << endl;

	delete[] times;
	delete[] solution;
}
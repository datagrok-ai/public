// testODE32sInGeneralSolverJNJ.cpp


#include<iostream>
#include<iomanip>
#include<ctime>
#include<fstream>
using namespace std;

#include <chrono>
using namespace std::chrono;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "odeSolver.h"
using namespace ode;

#include "tests.h"

namespace testODE32sGSjnj
{
	VectorXd F(double _time, VectorXd& y)
	{
		VectorXd res(y.size());

		// JNJ formulas are removed!
		return res;
	}

	
	MatrixXd J(double t, VectorXd& y, double eps)
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

	VectorXd T(double t, VectorXd& y, double eps)
	{
		return (F(t + eps, y) - F(t, y)) / eps;
	}


	VectorXd exactSolution(double t, unsigned dim)
	{
	// JNJ formulas are removed!

		return res;
	}
} // testODE32sGSjnj

// test of ode32s solver in generalized solver: multi-dimensional case, jnj problem
void testODE32sInGeneralSolverJNJ()
{
	using namespace testODE32sGSjnj;

	double t0 = 0.0;
	double t1 = 100.0;
		//1e+3;
	double h =
		0.00005;
		//1e-5;
		//0.0002;
	    // 1e-4;
	unsigned n = static_cast<unsigned>((t1 - t0) / h) + 1;
	unsigned d = 13;
	           // = 7;

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

	//auto start = time(0);
	auto start = high_resolution_clock::now();

	cout << "\nResult code: "
		<< fixedStepSolver(F, T, J, times, n, y.data(), d, solution)
		//<< oneStepSolver(F, times, n, y.data(), d, getNextPointRK4, solution)
		<< endl;

	//auto finish = time(0);
	auto finish = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(finish - start);

	cout << "\nTime is " << 1e-6 * duration.count() << " sec.\n";

	double mad = 0.0;

	VectorXd left(d);
	VectorXd right(d);
	VectorXd yPrev(d);
	VectorXd yNext(d);
	VectorXd yExact(d);
	
	for (unsigned i = 1; i < n - 1; i++)
	{
		// skip the point t = 120, since it's a discontinuity
		if ((times[i + 1] > 120) && (times[i - 1] < 120))
			continue;

		for (unsigned j = 0; j < d; j++)
			y(j) = solution[i + j * n];

		for (unsigned j = 0; j < d; j++)
			yPrev(j) = solution[i - 1 + j * n];

		for (unsigned j = 0; j < d; j++)
			yNext(j) = solution[i + 1 + j * n];

		right = F(times[i], y);

		left = (yNext - yPrev) / (times[i + 1] - times[i - 1]);		

		mad = fmax(mad, (left - right).cwiseAbs().maxCoeff());		

		//cout << "\n  iter no. " << i << "  " << (left - right).cwiseAbs().maxCoeff() << endl;
	}

	cout << "\nMAX ERROR: " << mad << endl;

	delete[] times;
	delete[] solution;
}
// testAdaptiveStepRKCK.cpp

#include<iostream>
#include<iomanip>
#include<ctime>
#include<fstream>
#include<vector>
#include<list>
#include <chrono>
using namespace std;
using namespace std::chrono;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "odeSolver.h"
using namespace ode;

#include "tests.h"

namespace testRKCK
{
	const unsigned DIM = 2;

	/*VectorXd f(double t, VectorXd& y)
	{
		VectorXd res(y.size());

		res(0) = 2.0 * t;
		res(1) = 3.0 * t * t;

		return res;
	}

	VectorXd exact(double t, unsigned dim)
	{
		VectorXd res(dim);

		res(0) = t * t;
		res(1) = t * t * t;

		return res;
	}*/

	//VectorXd f(double t, VectorXd& y)
	//{
	//	VectorXd res(y.size());

	//	res(0) = sin(t);
	//	res(1) = cos(t);

	//	return res;
	//}

	//VectorXd exact(double t, unsigned dim)
	//{
	//	VectorXd res(dim);

	//	res(0) = -cos(t);
	//	res(1) = sin(t);

	//	return res;
	//}

	//VectorXd f(double t, VectorXd& y)
	//{
	//	VectorXd res(y.size());

	//	res(0) = 2.0 * y(1);
	//	res(1) = 2.0 * y(0);

	//	return res;
	//}

	//VectorXd exact(double t, unsigned dim)
	//{
	//	VectorXd res(dim);

	//	res(0) = exp(2.0 * t) + exp(-2.0 * t);
	//	res(1) = exp(2.0 * t) - exp(-2.0 * t);

	//	return res;
	//}

	VectorXd f(double t, VectorXd& y)
	{
		VectorXd res(y.size());

		res(0) = -5.0 * y(0) + 3.0 * y(1);
		res(1) = 100.0 * y(0) - 301.0 * y(1);

		return res;
	}

	VectorXd exact(double t, unsigned dim)
	{
		VectorXd res(dim);

		res(0) = 52.96 * exp(-3.9899 * t) - 0.67 * exp(-302.0101 * t);
		res(1) = 17.83 * exp(-3.9899 * t) + 65.99 * exp(-302.0101 * t);

		return res;
	}
}; // testRKCK

namespace testRKCKjnj {

	const unsigned DIM = 13;

	VectorXd f(double _time, VectorXd& y)
	{
		VectorXd res(y.size());

		// JNJ formulas are removed!
		return res;
	}

	VectorXd exact(double t, unsigned dim)
	{
		// JNJ formulas are removed!

		return res;
	}
}; // testRKCKjnj

//using namespace testRKCK;
using namespace testRKCKjnj;


// test adaptive step Runge-Kutta (Cash-Karp) method 
void testAdaptiveStepRKCK()
{
	cout << "Test of generalized adaptive step R.-K. Cash-Karp solver.\n";

	unsigned dim = DIM;
	double t0 = 0.0;
	double t1 = 1000.0;	
	double hInitial = 0.01;
	double tol = 5e-5;
	VectorXd y0 = exact(t0, dim);

	/*list<double> times;
	list<VectorXd> solutions;*/

	vector<double> times;
	vector<VectorXd> solutions;

	cout << "\nSolvinig ...\n";

	auto start = high_resolution_clock::now();

	int resultCode = adaptiveStepSolver(f, t0, t1, hInitial, y0, tol, times, solutions);
	if (resultCode != NO_ERRORS)
	{
		cout << "\nFAIL! Result code: " << resultCode << endl;
		return;
	}

	cout << "\nSUCCESS!\n";
	
	auto finish = high_resolution_clock::now();


	//cout << "\n           t      y1(approx)      y2(approx)       y1(exact)       y2(exact)           error\n";

	//auto tIter = times.begin();
	//auto yIter = solutions.begin();

	//for (; tIter != times.end(); ++tIter, ++yIter)
	//{
	//	double t = *tIter;
	//	auto y = *yIter;
	//	auto yExact = exact(t, dim);
	//	double error = (y - yExact).cwiseAbs().maxCoeff();
	//	cout << setiosflags(ios::right)
	//		<< setw(12) << t
	//		<< setw(16) << y(0)
	//		<< setw(16) << y(1)
	//		<< setw(16) << yExact(0)
	//		<< setw(16) << yExact(1)
	//		<< setw(16) << error
	//		<< endl;
	//} // for

	auto duration = duration_cast<microseconds>(finish - start);

	cout << "\nTime of ODE solving is " << 1e-6 * duration.count() << " sec.\n";
}
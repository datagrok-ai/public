// testAdaptiveStepRKCK.cpp

#include<iostream>
#include<iomanip>
#include<ctime>
#include<fstream>
#include<vector>
#include<list>
#include <chrono>
#include<cstring>
#include<cmath>
using namespace std;
using namespace std::chrono;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

//#include "odeSolver.h"
//using namespace ode;

#include "odes.h"
using namespace odes;

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

	VectorXd f(double t, VectorXd& y)
	{
		VectorXd res(y.size());

		res(0) = 2.0 * y(1);
		res(1) = 2.0 * y(0);

		return res;
	}

	VectorXd exact(double t, unsigned dim)
	{
		VectorXd res(dim);

		res(0) = exp(2.0 * t) + exp(-2.0 * t);
		res(1) = exp(2.0 * t) - exp(-2.0 * t);

		return res;
	}

	MatrixXd J(double t, VectorXd& y, double eps)
	{
		MatrixXd res(y.size(), y.size());

		res << 0.0, 2.0,
			2.0, 0.0;

		return res;
	}

	VectorXd T(double t, VectorXd& y, double eps)
	{
		VectorXd res(y.size());

		res << 0.0,
			0.0;

		return y;
	}	
	
	/*
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
	}*/
}; // testRKCK

namespace testRKCKjnj {

	const unsigned DIM = 13;

	float yInitial[] = { // ORIGINAL INITIAL VALUES ARE REMOVED!
	};

	

		// Version 3
	VectorXd f(double _time, VectorXd& y) noexcept
	{
		VectorXd res(y.size());

		// ORIGINAL EXPRESSIONS ARE REMOVED
		
		return res;
	} // f

	MatrixXd J(double t, VectorXd& y, double eps)
	{
		MatrixXd res(y.size(), y.size());

		VectorXd val = f(t, y);

		VectorXd yDer = y;

		for (int i = 0; i < y.size(); i++)
		{
			yDer(i) += eps;

			res.col(i) = (f(t, yDer) - val) / eps;

			yDer(i) -= eps;
		}

		return res;
	}

	VectorXd T(double t, VectorXd& y, double eps)
	{
		return (f(t + eps, y) - f(t, y)) / eps;
	}


}; // testRKCKjnj

// test adaptive step Runge-Kutta (Cash-Karp) method 
void testAdaptiveStepRKCK()
{
	using namespace testRKCKjnj;

	cout << "Test of adaptive step R.-K. Cash-Karp solver.\n";

    int dim = DIM;
	float t0 = 0.0f;
	float t1 = 1000.0f;// 1000;
	float h = 0.001f;
	float tol = 0.00005f;
	
	int colCount = dim + 1;
	int rowCount = static_cast<int>((t1 - t0) / h + 1);
	float* dataFrame = new float[colCount * rowCount];

	cout << "\nSolving ...\n";

	auto start = high_resolution_clock::now();

	// basic version of the solver
	//int resultCode = solveODE(f, t0, t1, h, yInitial, tol, dataFrame, rowCount, colCount);

	// light version of the solver
	int resultCode = solveODE(f, T, J, t0, t1, h, yInitial, tol, dataFrame, rowCount, colCount);
	//int resultCode = solveODE(f, t0, t1, h, yInitial, tol, dataFrame, rowCount, colCount);

	auto finish = high_resolution_clock::now();

	auto duration = duration_cast<microseconds>(finish - start);

	cout << "\nTime of ODE solving is " << 1e-6 * duration.count() << " sec.\n";

	if (resultCode != NO_ERRORS)
	{
		cout << "\nFAILED! Result code: " << resultCode << endl;
		return;
	}

	cout << "\nSUCCESS!\n";

	string nameOfFile = "fae_" + to_string(t0) + "_" + to_string(t1) + "_" + to_string(h) + "_LIGHT_bdf.csv";
		
	ofstream resFile(nameOfFile, ios::out);
	if (!resFile)
	{
		cout << "\nFile creation fail!\n";
		return;
	}

	cout << "\nSaving to file ...\n";

	resFile << "t, FFox, KKox, FFred, KKred, Ffree, Kfree, FKred, FKox, MEAthiol_t, CO2, yO2P, Cystamine, VL\n";
	
	for (int i = 0; i < rowCount; i++)
	{
		resFile << dataFrame[i];

		for (int j = 1; j < colCount; j++)
			resFile << "," << dataFrame[i + j * rowCount];

		resFile << endl;
	}

	cout << "\nDone!\n";

	delete[] dataFrame;
} // testAdaptiveStepRKCK

void testOfSolverWithReducedComplexity()
{
	using namespace testRKCK;

	float t0 = 0.0f;
	float t1 = 10.0f;
	float h = 0.01f;
	float tol = 0.00005f;

	int colCount = DIM + 1;
	int rowCount = static_cast<int>((t1 - t0) / h + 1);
	float* dataFrame = new float[colCount * rowCount];

	double t = 0.0;
	auto y = exact(t, DIM);

	float* yInitial = new float[DIM];

	for (int i = 0; i < DIM; i++)
		yInitial[i] = y(i);

	cout << "\nTest of solver with reduced spatial complexity.\n";

	cout << "\nSolving ...\n";

	auto start = high_resolution_clock::now();

	int resultCode = solveODE(f, T, J, t0, t1, h, yInitial, tol, dataFrame, rowCount, colCount);

	auto finish = high_resolution_clock::now();

	auto duration = duration_cast<microseconds>(finish - start);

	cout << "\nTime of ODE solving is " << 1e-6 * duration.count() << " sec.\n";

	cout << "\n          time          y1_approx           y1_exact          y2_approx           y2_exact\n";

	for (int i = 0; i < rowCount; i++)
	{
		y = exact(t, DIM);

		cout << setiosflags(ios::right)
			<< setw(14) << t
			<< setw(19) << dataFrame[i + rowCount]			
			<< setw(19) << y(0)
			<< setw(19) << dataFrame[i + 2 * rowCount]
			<< setw(19) << y(1)
			<< endl;;

		t += h;
	}

	delete[] dataFrame;
	delete[] yInitial;
} // testOfSolverWithReducedComplexity
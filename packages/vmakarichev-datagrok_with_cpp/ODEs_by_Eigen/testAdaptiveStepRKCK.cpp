// testAdaptiveStepRKCK.cpp

#include<iostream>
#include<iomanip>
#include<ctime>
#include<fstream>
#include<vector>
#include<list>
#include <chrono>
#include<cstring>
using namespace std;
using namespace std::chrono;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "odeSolver.h"
using namespace ode;

#include "tests.h"

namespace testRKCKjnj {

	const unsigned DIM = 13;

	float yInitial[] = { 
		// initial conditions are removed
	};

	VectorXd f(double _time, VectorXd& y)
	{
		VectorXd res(y.size());

		// expressions are removed

		return res;
	}

}; // testRKCKjnj

//using namespace testRKCK;
using namespace testRKCKjnj;

// test adaptive step Runge-Kutta (Cash-Karp) method 
void testAdaptiveStepRKCK()
{
	cout << "Test of generalized adaptive step R.-K. Cash-Karp solver.\n";

    int dim = DIM;
	float t0 = 0;
	float t1 = 200;
	float h = 0.01f;
	float tol = 0.00005f;

	int colCount = dim + 1;
	int rowCount = static_cast<int>((t1 - t0) / h + 1);
	float* dataFrame = new float[colCount * rowCount];

	cout << "\nSolving ...\n";

	auto start = high_resolution_clock::now();

	int resultCode = solveODE(f, t0, t1, h, yInitial, tol, dataFrame, rowCount, colCount);

	auto finish = high_resolution_clock::now();

	auto duration = duration_cast<microseconds>(finish - start);

	cout << "\nTime of ODE solving is " << 1e-6 * duration.count() << " sec.\n";

	if (resultCode != NO_ERRORS)
	{
		cout << "\nFAILED! Result code: " << resultCode << endl;
		return;
	}

	cout << "\nSUCCESS!\n";

	string nameOfFile = "fae_" + to_string(t0) + "_" + to_string(t1) + "_" + to_string(h) + ".csv";
	
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
}
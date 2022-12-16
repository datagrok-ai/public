// test.cpp
// Implementations of functions for testing Solver

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
using namespace std;

#include "test.h"

#include "solver.h"

#include "wrapper.h"

// try solver function
void trySolver()
{
	double t0 = 0;
	double t1 = 100;
	int N = 1000;
	double step = (t1 - t0) / N;

	vector<vector<double>> result;

	vector<double> times(N + 1);

	times[0] = t0;

	for (int i = 1; i <= N; i++)
		times[i] = times[i - 1] + step;

	/*cout << "\ntimes:\n";
	for (int i = 0; i <= N; i++)
	{
		cout << "  " << i << "  " << times[i] << endl;
	}*/

	vector<double> _parameters(1);
	_parameters[0] = 6.6;

	Solve(&result, &times, &_parameters, NULL, NULL, NULL);

	double VLinitial = _parameters[0];	
	double k1red = 0.000934 * 60;
	double k1ox = 0.00018 * 60;
	double MEAthiol = 34;
	double pH = 7.4;
	double pKa2MEA = 8.19;
	double MEAthiolate = MEAthiol * pow(10, pH - pKa2MEA);

	double c = -k1red * pow(VLinitial / 6.6, 3) * pow(MEAthiolate, 2);

	cout << "\n  time       FFred       KKred          VL        FFox        KKox        Analytic\n";

	for (int i = 0; i <= N; i++)
		cout << setiosflags(ios::right) << setw(6) << times[i] 
		     << setw(12) << result[0][i] 
		     << setw(12) << result[1][i] 
		     << setw(12) << result[4][i] 
		     << setw(12) << result[2][i]
		     << setw(12) << result[3][i]
		     << setw(16) << 50 * exp(c * times[i])
		     << endl;
}

// test of solver wrapper
void testSolverWrapper()
{
	float t0 = 0;
	float t1 = 100;
	int N = 1000;
	float step = (t1 - t0) / N;

	int tValsLength = N + 1;
	float * tVals = new float[tValsLength];

	tVals[0] = t0;
	for (int i = 1; i < tValsLength; i++)
		tVals[i] = tVals[i - 1] + step;

	/*cout << "Time:\n";
	for (int i = 0; i < tValsLength; i++)
		cout << "  " << tVals[i] << endl;*/

	int paramsLength = 1;
	float * params = new float[paramsLength];

	params[0] = 6.6f;

	int solutionRowCount = tValsLength;
	int solutionColumnCount = 5;

	float * solution = new float[tValsLength * 5];

	cout << "\nCall wrapper result: "
		<< solverWrapper(tVals, tValsLength, params, paramsLength, solution, solutionRowCount, solutionColumnCount)
		<< endl;

	cout << "\n  time       FFred       KKred          VL        FFox        KKox\n";

	for (int i = 0; i <= N; i++)
		cout << setiosflags(ios::right) << setw(6) << tVals[i]
		<< setw(12) << solution[i]
		<< setw(12) << solution[i + solutionRowCount]
		<< setw(12) << solution[i + 4 * solutionRowCount]
		<< setw(12) << solution[i + 2 * solutionRowCount]
		<< setw(12) << solution[i + 3 * solutionRowCount]
		<< endl;

	delete[] tVals;
	delete[] params;
	delete[] solution;
}
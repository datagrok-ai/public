// analyzeStiffProblemMultiDim.cpp

// Analyzing solver for stiff problem: one-dimensional case

#include<iostream>
#include<iomanip>
#include<ctime>
#include<fstream>
using namespace std;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "odeSolver.h"
using namespace ode;

// STIFF EXAMPLE: Steve Chapra and Raymond P. Canale. Numerical Methods for Engineers, 2021 (page 770).

VectorXd w(double t, VectorXd & y)
{
	VectorXd res(y.size());

	res(0) = -5.0 * y(0) + 3.0 * y(1);
	res(1) = 100.0 * y(0) - 301.0 * y(1);

	return res;
}

VectorXd exactSolut(double t, unsigned dim)
{
	VectorXd res(dim);

	res(0) = 52.96 * exp(-3.9899 * t) - 0.67 * exp(-302.0101 * t);
	res(1) = 17.83 * exp(-3.9899 * t) + 65.99 * exp(-302.0101 * t);

	return res;
}

// analysis of solving stiff problem: one-dimensional case
void analyzeStiffProblemMultiDim()
{
	cout << "Analysis of solver for stiff problem: multi-dimensional case.\n";

	unsigned dim = 2;

	double h = 0.0005;

	double t0 = 0.0;
	double t1 = 100.0;

	unsigned N = static_cast<unsigned>((t1 - t0) / h) + 1;

	cout << "\nThe initial problem is solved on the segment [t0, t1] with the step h,\nwhere\n"
		<< "   t0 = " << t0
		<< "\n   t1 = " << t1
		<< "\n   h = " << h
		<< "\n\n   number of points = " << N
		<< endl;

	double * times = new double[N];
	double * solution = new double[N * dim];

	times[0] = 0.0;
	VectorXd yInitial = exactSolut(times[0], dim);

	cout << "\nInitialization of the array 'times' ...\n";
	for (unsigned i = 1; i < N; i++)
		times[i] = times[i - 1] + h;

	cout << "\nSolving ODE ...\n";

	time_t start = time(0);

	int resCode = oneStepSolver(w, times, N, yInitial.data(), dim, getNextPointRK4, solution);

	time_t finish = time(0);

	cout << ((resCode == NO_ERRORS) ? "\nSUCCESS!" : "\nFAULT!") << endl;

	cout << "\nComputation time " << finish - start << " sec.\n";

	if (resCode != NO_ERRORS)
	{
		cout << "\nResult code: " << resCode << endl;
		return;
	}

	cout << "\nComputing error ...\n";

	double mae = fabs(solution[0] - yInitial(0));
	mae = fmax(mae, fabs(solution[N] - yInitial(1)));

	ofstream file("RK4.csv", ios::out);

	for (unsigned i = 1; i < N; i++)
	{
		yInitial = exactSolut(times[i], dim);
		mae = fmax(mae, fabs(solution[i] - yInitial(0)));
		mae = fmax(mae, fabs(solution[i + N] - yInitial(1)));
		file << mae << "\n";
	}

	cout << "\nMaximum absolute error: " << mae << endl;

	/*cout << "\n   Time   Solution(approx)   Solution(exact)   Error(absolute)   Error(relative)\n";
	for (unsigned i = 0; i < N; i++)
	cout << setw(7) << setiosflags(ios::right) << times[i]
	<< setw(19) << solution[i]
	<< setw(18) << exactSol(times[i])
	<< setw(18) << fabs(solution[i] - exactSol(times[i]))
	<< setw(18) << fabs((solution[i] - exactSol(times[i])) / exactSol(times[i]))
	<< endl;*/

	delete[] times;
	delete[] solution;

}
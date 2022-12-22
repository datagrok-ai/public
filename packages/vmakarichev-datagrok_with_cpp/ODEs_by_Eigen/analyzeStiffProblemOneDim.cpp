// analyzeStiffProblemOneDim.cpp

// Analyzing solver for stiff problem: one-dimensional case

#include<iostream>
#include<iomanip>
#include<ctime>
using namespace std;

#include "odeSolver.h"
using namespace ode;

// STIFF EXAMPLE: Steve Chapra and Raymond P. Canale. Numerical Methods for Engineers, 2021 (page 767).

double g(double t, double & y)
{
	return -1000.0 * y + 3000.0 - 2000.0 * exp(-t);
}

double exactSol(double t)
{
	return 3 - 0.998 * exp(-1000.0 * t) - 2.002 * exp(-t);
}

// analysis of solving stiff problem: one-dimensional case
void analyzeStiffProblemOneDim()
{
	cout << "Analysis of solver for stiff problem: one-dimensional case.\n";
	
	double h = 0.0002;

	double t0 = 0.0;
	double t1 = 10000.0;

	unsigned N = static_cast<unsigned>((t1 - t0) / h) + 1;

	cout << "\nThe initial problem is solved on the segment [t0, t1] with the step h,\nwhere\n"
		<< "   t0 = " << t0
		<< "\n   t1 = " << t1
		<< "\n   h = " << h
		<< "\n\n   number of points = " << N
		<< endl;

	double * times = new double[N];
	double * solution = new double[N];

	times[0] = 0.0;
	double yInitial = exactSol(times[0]);

	cout << "\nInitialization of the array 'times' ...\n";
	for (unsigned i = 1; i < N; i++)
		times[i] = times[i - 1] + h;

	cout << "\nSolving ODE ...\n";

	time_t start = time(0);

	int resCode = oneStepSolver(g, times, N, yInitial, getNextPointRK5, solution);

	time_t finish = time(0);

	cout << ( (resCode == NO_ERRORS) ? "\nSUCCESS!" : "\nFAULT!" )<< endl;		

	cout << "\nComputation time " << finish - start << " sec.\n";

	if (resCode != NO_ERRORS)
	{
		cout << "\nResult code: " << resCode << endl;
		return;
	}

	cout << "\nComputing error ...\n";

	double mae = fabs(solution[0] - exactSol(times[0]));

	for (unsigned i = 0; i < N; i++)
	{
		mae = fmax(mae, fabs(solution[i] - exactSol(times[i])));
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
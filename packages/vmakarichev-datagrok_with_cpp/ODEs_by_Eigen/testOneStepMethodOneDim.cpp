// testOneStepMethodOneDim.cpp

// Implementation of the function testOneStepMethodOneDim() - test of one-step method: one-dimensional case

#include<iostream>
#include<iomanip>
using namespace std;

#include "odeSolver.h"
using namespace ode;

/*
double f(double t, double & y)
{
	return 2 * t;
}

double exactSolution(double t)
{
	return t * t;
}*/

/*double f(double t, double & y)
{
	return cos(t);
}

double exactSolution(double t)
{
	return sin(t);
}*/

/*double f(double t, double & y)
{
	return exp(t);
}

double exactSolution(double t)
{
	return exp(t);
}*/

//double f(double t, double & y)
//{
//	return 2.0 * y - 2.0;
//}
//
//double exactSolution(double t)
//{
//	return exp(2.0 * t) + 1.0;
//}

//double f(double t, double & y)
//{
//	return 10.0 * y - 2.0;
//}
//
//double exactSolution(double t)
//{
//	return exp(10.0 * t) + 0.2;
//}

//double f(double t, double & y)
//{
//	return -10.0 * y - 2.0;
//}
//
//double exactSolution(double t)
//{
//	return exp(-10.0 * t) - 0.2;
//}

// STIFF EXAMPLE: Steve Chapra and Raymond P. Canale. Numerical Methods for Engineers, 2021 (page 767).

double f(double t, double & y)
{
	return -1000.0 * y + 3000.0 - 2000.0 * exp(-t);
}

double exactSolution(double t)
{
	return 3 - 0.998 * exp(-1000.0 * t) - 2.002 * exp(-t);
}

// test of one-step method: one-dimensional case
void testOneStepMethodOneDim()
{
	cout << "Test of one step method: one-dimensional case.\n";

	/*double t = 0.0;
	double h = 0.001;
	double y = 0.0;

	int n = 100;

	cout << "\nSolution:\n";
	cout << "  i = " << 0 << "   t = " << t << "   y = " << y << endl;

	for (int i = 1; i <= n; i++)
	{
		y = getNextPointRK4(f, t, y, h);
		t += h;

		cout << "  i = " << i << "   t = " << t << "   y = " << y << endl;
	}*/

	double h = 0.0005;	
	const int N = 10000001;

	double * times= new double [N];
	double * solution = new double [N];

	times[0] = 0.0;
	double yInitial = exactSolution(times[0]);

	for (int i = 1; i < N; i++)
		times[i] = times[i - 1] + h;

	cout << "Solving ODE result code: "
		<< oneStepSolver(f, times, N, yInitial, getNextPointRK4, solution)  << endl;

	cout << "\n   Time   Solution(approx)   Solution(exact)   Error(absolute)   Error(relative)\n";
	for (int i = 0; i < N; i++)
		cout << setw(7) << setiosflags(ios::right) << times[i]
		<< setw(19) << solution[i]
		<< setw(18) << exactSolution(times[i])
		<< setw(18) << fabs(solution[i] - exactSolution(times[i]))
		<< setw(18) << fabs((solution[i] - exactSolution(times[i])) / exactSolution(times[i]))
		<< endl;

	delete[] times;
	delete[] solution;

} // testRK4oneDim
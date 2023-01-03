// testOfStiffExampleDemo.cpp

#include<iostream>
#include<iomanip>
#include<ctime>
#include<fstream>
using namespace std;

#include "tests.h"

#include "demos.h"

void testOfStiffExampleDemo()
{
	int dim = 2;
	float t0 = 0.0;
	float t1 = 1.0;
	float step = 0.01;

	int resultRowCount = static_cast<int>((t1 - t0) / step) + 1;
	int resultColCount = dim + 1;
	float* result = new float[resultRowCount * resultColCount];

	cout << "\nSolving ...\n";

	cout << "\nResult code: "
		<< stiffExample(t0, t1, step, result, resultRowCount, resultColCount)
		<< endl;

	cout << "    time        y(1)        y(2)\n";
	for (int i = 0; i < resultRowCount; i++)
		cout << setiosflags(ios::right)
		<< setw(8) << result[i]
		<< setw(12) << result[i + resultRowCount]
		<< setw(12) << result[i + 2 * resultRowCount]
		<< endl;

	delete[] result;
}

void testOfStiffExampleRKDemo()
{
	int dim = 2;
	float t0 = 0.0;
	float t1 = 1.0;
	float step = 0.001;

	int resultRowCount = static_cast<int>((t1 - t0) / step) + 1;
	int resultColCount = dim + 1;
	float* result = new float[resultRowCount * resultColCount];

	cout << "\nSolving ...\n";

	cout << "\nResult code: "
		<< stiffExampleRK(t0, t1, step, result, resultRowCount, resultColCount)
		<< endl;

	cout << "    time        y(1)        y(2)\n";
	for (int i = 0; i < resultRowCount; i++)
		cout << setiosflags(ios::right)
		<< setw(8) << result[i]
		<< setw(12) << result[i + resultRowCount]
		<< setw(12) << result[i + 2 * resultRowCount]
		<< endl;

	delete[] result;
}
#include<iostream>
using namespace std;

#include "test.h"
#include "fixedPoint.h"

vector<double> func(double t, vector<double> & y)
{
	vector<double> res(y.size());

	res[0] = y[0] + y[1];
	res[1] = y[0] - y[1];

	return res;
}

// test computation of k-s using fixedPoint
void testFixedPoint()
{
	int N = 2;
	int s = 4;

	double t = 1.3;

	vector<double> y(N);
	y[0] = 2.1;
	y[1] = 1.2;

	vector<vector<double>> k(N);
	for (int i = 0; i < N; i++)
		k[i].resize(s);

	cout << "\nk:\n";
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < s; j++)
			cout << "  " << k[i][j];
		cout << endl;
	}

	vector<vector<double>> a(s);
	for (int i = 0; i < s; i++)
		a[i].resize(s);

	a[0][0] = 0.1; a[0][1] = 0.2; a[0][2] = 0.3; a[0][3] = 0.2;
	a[1][0] = 0.2; a[1][1] = 0.3; a[1][2] = 0.2; a[1][3] = 0.1;
	a[2][0] = 0.1; a[2][1] = 0.1; a[2][2] = 0.2; a[2][3] = 0.1;
	a[3][0] = 0.1; a[3][1] = 0.2; a[3][2] = 0.1; a[3][3] = 0.1;

	cout << "\na:\n";
	for (int i = 0; i < s; i++)
	{
		for (int j = 0; j < s; j++)
			cout << "  " << a[i][j];
		cout << endl;
	}

	vector<double> c(s);
	for (int i = 0; i < s; i++)
	{
		c[i] = a[i][0];

		for (int j = 1; j <= i; j++)
			c[i] += a[i][j];
	}

	cout << "\nc:\n";
	for (int i = 0; i < s; i++)
		cout << "  " << c[i];
	cout << endl;

	double h = 0.1;

	double precision = 0.000001;

	getRoots(func, t, y, k, a, c, h, precision);

	cout << "\nFinal result k:\n";
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < s; j++)
			cout << "  " << k[i][j];
		cout << endl;
	}
}




// test computation of k-s using fixedPoint, another example
void testFixedPointAnother()
{
	int N = 7;
	int s = 1;

	double t = 2.331;

	vector<double> y(N);
	y[0] = 2.1;
	y[1] = 1.2;
	y[2] = 1.2;
	y[3] = 2.1;
	y[4] = 1.2;
	y[5] = 1.2;
	y[6] = 1.2;

	vector<vector<double>> k(N);
	for (int i = 0; i < N; i++)
		k[i].resize(s);

	cout << "\nk:\n";
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < s; j++)
			cout << "  " << k[i][j];
		cout << endl;
	}

	vector<vector<double>> a(s);
	for (int i = 0; i < s; i++)
		a[i].resize(s);

	a[0][0] = 1;

	cout << "\na:\n";
	for (int i = 0; i < s; i++)
	{
		for (int j = 0; j < s; j++)
			cout << "  " << a[i][j];
		cout << endl;
	}

	vector<double> c(s);
	for (int i = 0; i < s; i++)
	{
		c[i] = a[i][0];

		for (int j = 1; j <= i; j++)
			c[i] += a[i][j];
	}

	cout << "\nc:\n";
	for (int i = 0; i < s; i++)
		cout << "  " << c[i];
	cout << endl;

	double h = 0.00001;

	double precision = 0.01;

	getRoots(DerivFunc, t, y, k, a, c, h, precision);

	cout << "\nFinal result k:\n";
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < s; j++)
			cout << "  " << k[i][j];
		cout << endl;
	}
}

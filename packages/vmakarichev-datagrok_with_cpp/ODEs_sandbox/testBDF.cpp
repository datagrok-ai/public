// testBDF.cpp

#include<iostream>
#include<vector>
using namespace std;

typedef vector<double> Vector;

#include "test.h"
#include "BDF.h"

Vector f(double t, Vector & y)
{
	Vector result(y.size());

	result[0] = y[0] * y[1] + y[1] * y[2] + y[0] * y[2];
	result[1] = y[0] * y[0] + y[1] * y[1] + y[2] * y[2];
	result[2] = y[0] * y[0] * y[0] + y[1] * y[1] * y[1] + y[2] * y[2] * y[2];

	/*result[0] = y[0] + y[1] + y[2];
	result[1] = 2 * y[0] - y[1] + 3 * y[2];
	result[2] = y[0] - 2 * y[1] + 2 * y[2];*/

	return result;
}

// test of BDF: one step computation
void testBDF()
{
	double t = 1.1;
	double h = 0.0031;

	const unsigned BDF_MAX_STEP = 1;
	int N = 3;

	vector<Vector *> yFoundPtrs;
	yFoundPtrs.reserve(BDF_MAX_STEP);

	cout << yFoundPtrs.capacity() << endl << yFoundPtrs.size() << endl;

	Vector y0(N);
	y0[0] = 1.1;
	y0[1] = -0.34;
	y0[2] = 0.311;

	/*y0[0] = 1.1; 
	y0[1] = -2.34;
	y0[2] = 31.1;*/

	Vector y1(N);
	y1[0] = 11.13;
	y1[1] = -0.34;
	y1[2] = 21.1;

	Vector y2(N);
	y2[0] = 7.11;
	y2[1] = 2.134;
	y2[2] = 11.1;

	Vector y3(N);
	y3[0] = 2.1;
	y3[1] = 4.134;
	y3[2] = 0.1;

	Vector y4(N);
	y4[0] = -2.31;
	y4[1] = 11.134;
	y4[2] = -5.51;

	Vector y5(N);
	y5[0] = -5.31;
	y5[1] = 19.134;
	y5[2] = -15.51;

	yFoundPtrs.push_back(&y0);
	yFoundPtrs.push_back(&y1);
	yFoundPtrs.push_back(&y2);
	yFoundPtrs.push_back(&y3);
	yFoundPtrs.push_back(&y4);
	yFoundPtrs.push_back(&y5);

	cout << yFoundPtrs.capacity() << endl << yFoundPtrs.size() << endl;

	auto res = f(t, y0);

	cout << "\nres:\n";
	for (int i = 0; i < N; i++)
		cout << "  " << res[i];

	auto y = bdf(f, t, h, yFoundPtrs);

	cout << "\ny:\n";
	for (int i = 0; i < N; i++)
		cout << "  " << y[i];

} // testBDF
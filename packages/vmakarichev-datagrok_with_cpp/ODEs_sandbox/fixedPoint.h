// fixedPoint.h

#ifndef FIXED_POINT_H
#define FIXED_POINT_H

#include<vector>
using std::vector;

/* Computation of the coefficients k_{i,j}, required for non-explicit Runge-Kutta method.
func - right part of the system of ODEs solved, i.e. dy/dt = func(t,y);
t - value of t[n], previous time;
y - values of y(t[n]) <- these values are computed at the previous step of R.-K. method;
k - coefficients { k_{i,j} } that are computed;
a, c - parameters of R.-K. method;
h - step (h = t[n + 1] - t[n]);
precision - precision of the applied iteration approach.
*/
void getRoots(vector<double>(*func) (double, vector<double> &),
	double t,
	vector<double> & y,
	vector<vector<double>> & k,
	vector<vector<double>> & a,
	vector<double> & c,
	double h,
	double precision);

#endif // !FIXED_POINT_H



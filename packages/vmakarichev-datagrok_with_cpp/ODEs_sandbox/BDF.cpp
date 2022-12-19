// BDF.cpp

#include<iostream>
using namespace std;

#include "BDF.h"

/* Returns solution of the system y = c * h * f(t, y) + v, 
   where the number c and the vector v are defined by BDF formulas (see https://en.wikipedia.org/wiki/Backward_differentiation_formula).
      f - vector-function;
	  t - time;
	  h - step;
	  yFoundPtrs - pointers to solution of the OEDs system at the previos steps;
	  precision - precision of solution;
	  eps - step, which is used in order to estimate partial derivatives.
*/
std::vector<double> bdf(std::vector<double>(*f) (double, std::vector<double> &),
	double t,
	double h,
	std::vector<std::vector<double> *> & yFoundPtrs,
	double precision,
	double eps)
{
	// step of BDF (see formulas)
	unsigned s = yFoundPtrs.size(); 

	// dimension of vectors, which are operated; also, a number of ODEs
	unsigned N = yFoundPtrs[0]->size(); 

	cout << "\nN = " << N << endl;
	cout << "\ns = " << s << endl;
	
	// constant vector, defined by BDF step
	std::vector<double> v(N);

	/*// coumputation of the constant c and the vector v (see https://en.wikipedia.org/wiki/Backward_differentiation_formula)
	switch (s)
	{
	case 1: // BSD1
		c = 1.0; 

		for (unsigned i = 0; i < N; i++)
			v[i] = -yFoundPtrs[0]->at(i);

		cout << "\nv:\n";
		for (unsigned i = 0; i < N; i++)
			cout << "  " << v[i];

		break;

	case 2: // BSD2
		c = 2.0 / 3.0;

		for (unsigned i = 0; i < N; i++)
			v[i] = -yFoundPtrs[1]->at(i) * 4.0 / 3.0 + yFoundPtrs[0]->at(i) / 3.0;

		break;

	case 3: // BSD3
		c = 6.0 / 11.0;

		for (unsigned i = 0; i < N; i++)
			v[i] = -yFoundPtrs[2]->at(i) * 18.0 / 11.0 + yFoundPtrs[1]->at(i) * 9.0 / 11.0 - yFoundPtrs[0]->at(i) * 2.0 / 11.0;

		break;

	case 4: // BSD4
		c = 12.0 / 25.0;

		for (unsigned i = 0; i < N; i++)
			v[i] = -yFoundPtrs[3]->at(i) * 48.0 / 25.0 + yFoundPtrs[2]->at(i) * 36.0 / 25.0
			- yFoundPtrs[1]->at(i) * 16.0 / 25.0 + yFoundPtrs[0]->at(i) * 3.0 / 25.0;

		break;

	default:
		break; // TODO: add processing this case; actually, this is a mistake!
	}*/

	double prod = _c[s - 1] *  h;

	for (unsigned i = 0; i < N; i++)
	{
		v[i] = _mult[s - 1][0] * yFoundPtrs[0]->at(i);

		for (unsigned j = 1; j < s; j++)
			v[i] += _mult[s - 1][j] * yFoundPtrs[j]->at(i);
	}

	// initialization (the last of y's, which has been already found, is used)
	std::vector<double> yPrev(* yFoundPtrs[s - 1]);

	cout << "\nyPrev:\n";
	for (unsigned i = 0; i < N; i++)
		cout << "  " << yPrev[i];

	// solution of the problem considered
	std::vector<double> y(N);

	// maximum absolute deviation
	double mad = 0.0;

	int iter = 1;

	// itearative search for the solution
	while (true)
	{
		std::vector<double> res = f(t, yPrev);

		// computing coordinates of y
		for (unsigned i = 0; i < N; i++)
		{
			std::vector<double> yDer(yPrev);
			yDer[i] += eps;

			std::vector<double> resDer = f(t, yDer);

			double lambda = (prod * (res[i] - resDer[i]) - eps) / eps;

			// TODO: add processing the case lambda = 0

			y[i] = yPrev[i] - (prod * res[i] - v[i] - yPrev[i]) / lambda;
		} // for i

		// check precision
		mad = fabs(y[0] - yPrev[0]);

		for (unsigned i = 1; i < N; i++)
			mad = fmax(mad, fabs(y[i] - yPrev[i]));

		cout << "\n iter no. " << iter++ << "   mad = " << mad << endl;

		if (mad < precision)
			break;
		else
			for (unsigned i = 0; i < N; i++)
				yPrev[i] = y[i];
	} // while

	// check result: the difference c*h*f(t,y)-v-y MUST BE small. 
	std::vector<double> res = f(t, y);

	mad = fabs(prod * res[0] - v[0] - y[0]);
	for (unsigned i = 1; i < N; i++)
		mad = fmax(mad, fabs(prod * res[i] - v[i] - y[i]));

	cout << "\nError: mad = " << mad << endl;

	return y;
}// bdf
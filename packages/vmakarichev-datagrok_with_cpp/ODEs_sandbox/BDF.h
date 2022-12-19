// BDF.h

#ifndef BDF_H
#define BDF_H

#include<vector>

const unsigned MAX_BDF_STEP = 6;

// special quantities specified in https://en.wikipedia.org/wiki/Backward_differentiation_formula
const double _c[MAX_BDF_STEP] = {1.0, 2.0 / 3.0, 6.0 / 11.0, 12.0 / 25.0, 60.0 / 137.0, 60.0 / 147.0};

// special quantities specified in https://en.wikipedia.org/wiki/Backward_differentiation_formula
const double _mult[MAX_BDF_STEP][MAX_BDF_STEP] = { -1.0, 0, 0, 0, 0, 0,
                                                    1.0 / 3.0, -4.0 / 3.0, 0, 0, 0, 0,
											       -2.0 / 11.0, 9.0 / 11.0, -18.0 / 11.0, 0, 0, 0,
										      	    3.0 / 25.0, -16.0 / 25.0, 36.0 / 25.0, -48.0 / 25.0, 0, 0,
											      -12.0 / 137.0, 75.0 / 137.0, -200.0 / 137.0, 300.0 / 137.0, -300.0 / 137.0, 0,
											       10.0 / 147.0, -72.0 / 147.0, 225.0 / 147.0, -400.0 / 147.0, 450.0 / 147.0, -360.0 / 147.0 };

/* Returns solution of the system y = c * h * f(t, y) - v, 
   where the number c and the vector v are defined by BDF formulas (see https://en.wikipedia.org/wiki/Backward_differentiation_formula).
      f - vector-function;
      t - time;
      h - step;
      yFoundPtrs - pointers to solution of the OEDs system at the previos steps;
	  precision - precision of solution;
	  eps - step, which is used in order to estimate partial derivatives.
*/
std::vector<double> bdf(std::vector<double> (*f) (double, std::vector<double> &), 
	double t,
	double h,
	std::vector<std::vector<double> *> & yFoundPtrs,
	double precision = 0.001,
	double eps = 0.001);

#endif // !BDF_H


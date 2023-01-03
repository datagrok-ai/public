// tryAdaptiveStepRKCKmethod.cpp

#include<iostream>
#include<iomanip>
#include<ctime>
#include<fstream>
using namespace std;

#include <chrono>
using namespace std::chrono;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "tests.h"

namespace tryRKCK
{
	const unsigned DIM = 2;

	/*VectorXd f(double t, VectorXd& y)
	{
		VectorXd res(y.size());

		res(0) = 2.0 * t;
		res(1) = 3.0 * t * t;

		return res;
	}

	VectorXd exact(double t, unsigned dim)
	{
		VectorXd res(dim);

		res(0) = t * t;
		res(1) = t * t * t;

		return res;
	}*/

	//VectorXd f(double t, VectorXd& y)
	//{
	//	VectorXd res(y.size());

	//	res(0) = sin(t);
	//	res(1) = cos(t);

	//	return res;
	//}

	//VectorXd exact(double t, unsigned dim)
	//{
	//	VectorXd res(dim);

	//	res(0) = -cos(t);
	//	res(1) = sin(t);

	//	return res;
	//}

	//VectorXd f(double t, VectorXd& y)
	//{
	//	VectorXd res(y.size());

	//	res(0) = 2.0 * y(1);
	//	res(1) = 2.0 * y(0);

	//	return res;
	//}

	//VectorXd exact(double t, unsigned dim)
	//{
	//	VectorXd res(dim);

	//	res(0) = exp(2.0 * t) + exp(-2.0 * t);
	//	res(1) = exp(2.0 * t) - exp(-2.0 * t);

	//	return res;
	//}

	VectorXd f(double t, VectorXd& y)
	{
		VectorXd res(y.size());

		res(0) = -5.0 * y(0) + 3.0 * y(1);
		res(1) = 100.0 * y(0) - 301.0 * y(1);

		return res;
	}

	VectorXd exact(double t, unsigned dim)
	{
		VectorXd res(dim);

		res(0) = 52.96 * exp(-3.9899 * t) - 0.67 * exp(-302.0101 * t);
		res(1) = 17.83 * exp(-3.9899 * t) + 65.99 * exp(-302.0101 * t);

		return res;
	}
}; // tryRKCK

namespace tryRKCKjnj {

	const unsigned DIM = 13;

	VectorXd f(double _time, VectorXd& y)
	{
		VectorXd res(y.size());

		// JNJ formulas are removed!
		return res;
	}

	VectorXd exact(double t, unsigned dim)
	{
		VectorXd res(dim);
		// JNJ formulas are removed!

		return res;
	}
}; // tryRKCKjnj

//using namespace tryRKCK;
using namespace tryRKCKjnj;

int RKCK(VectorXd& y, VectorXd& dydt, double t, double h, VectorXd& yOut, VectorXd& yErr)
{
	VectorXd ytemp = y + 0.2 * h * dydt;
	VectorXd k2 = f(t + 0.2 * h, ytemp);
	ytemp = y + h * (0.075 * dydt + 0.225 * k2);
	VectorXd k3 = f(t + 0.3 * h, ytemp);
	ytemp = y + h * (0.3 * dydt - 0.9 * k2 + 1.2 * k3);
	VectorXd k4 = f(t + 0.6 * h, ytemp);
	ytemp = y + h * (-0.2037037037037037 * dydt + 2.5 * k2 - 2.592592592592593 * k3 + 1.296296296296296 * k4);
	VectorXd k5 = f(t + 0.6 * h, ytemp);
	ytemp = y + h * (0.0294958043981481 * dydt + 0.341796875 * k2 + 0.0415943287037037 * k3 
		+ 0.4003454137731481 * k4 + 0.061767578125 * k5);
	VectorXd k6 = f(t + 0.875 * h, ytemp);
	yOut = y + h * (0.0978835978835979 * dydt + 0.4025764895330113 * k3 
		+ 0.21043771043771045 * k4 + 0.2891022021456804 * k6);
	yErr = h * (-0.004293774801587311 * dydt + 0.018668586093857853 * k3 - 0.034155026830808066 * k4 
		- 0.019321986607142856 * k5 + 0.03910220214568039 * k6);

	return 0;
} // RKCK

int adapt(VectorXd& y, VectorXd& dydt, double & t, double hTry, VectorXd& yScale, double eps, double & hNext)
{
	const double SAFETY = 0.9;
	const double PSHRNK = -0.25;
	const double PSGROW = -0.2;
	const double REDUCE_COEF = 0.25;
	const double GROW_COEF = 4.0;
	const double ERR_CONTR = 1.89e-4;

	double h = hTry;	

	unsigned dim = y.size();

	VectorXd yTemp(dim);
	VectorXd yErr(dim);

	while (true)
	{
		RKCK(y, dydt, t, h, yTemp, yErr);

		double errmax = 0.0;

		for (unsigned i = 0; i < dim; i++)
			errmax = fmax(errmax, fabs(yErr(i) / yScale(i)));

		errmax /= eps;

		if (errmax > 1)
		{
			double hTemp = SAFETY * h * pow(errmax, PSHRNK);

			h = fmax(hTemp, REDUCE_COEF * h);

			double tNew = t + h;

			if (tNew == t)
				return -1;
		}
		else
		{
			if (errmax > ERR_CONTR)
				hNext = SAFETY * h * pow(errmax, PSGROW);
			else
				hNext = GROW_COEF * h;

			t = t + h;
			y = yTemp;
			break;
		}
	}

	return 0;
} // adapt

// try adaptive step Runge-Kutta (Cash-Karp) method 
void tryAdaptiveStepRKCKmethod()
{
	unsigned dim = DIM;
	double t0 = 0.0;
	double t1 = 10.0;
	VectorXd y0 = exact(t0, dim);

	double hInitial = 0.01;
	double eps = 5e-5;	
	VectorXd tiny = VectorXd::Ones(dim) * 1e-20;

	//cout << "\ntiny = " << tiny.transpose() << endl;

	double t = t0;
	double h = hInitial;
	double hNext = 0.0;
	VectorXd y = exact(t0, dim);
	VectorXd dydt(dim);
	VectorXd yScale(dim);
	VectorXd yPrev(dim);

	//cout << "\n t = " << t << "  y = " << y.transpose() << endl;

	unsigned iter = 1;

	bool flag = true;

	double error = 0.0;
	double errApprox = 0.0;

	//cout << "\n           t      y1(approx)      y2(approx)       y1(exact)       y2(exact)           error     errorApprox\n";

	auto start = high_resolution_clock::now();	

	double hMax = 0.0;

	while (flag)
	{
		dydt = f(t, y);

		yPrev = y;

		yScale = y.cwiseAbs() + (dydt * h).cwiseAbs() + tiny;

		if (t + h > t1) {
			h = t1 - t;
			flag = false;
		}

		/*cout << "\n  y' = " << dydt.transpose() << endl;

		cout << "\n |y| = " << y.cwiseAbs() << endl;

		cout << "\n yScale = " << yScale.transpose() << endl;*/
		
		adapt(y, dydt, t, h, yScale, eps, hNext);

		errApprox = fmax(errApprox, ((y - yPrev) / h - dydt).cwiseAbs().maxCoeff());

		h = hNext;

		hMax = fmax(h, hMax);

		/*VectorXd yExact = exact(t, dim);

		error = (y - yExact).cwiseAbs().maxCoeff();			

		cout << setiosflags(ios::right)
			<< setw(12) << t
			<< setw(16) << y(0)
			<< setw(16) << y(1)
			<< setw(16) << yExact(0)
			<< setw(16) << yExact(1)
			<< setw(16) << error
			<< setw(16) << errApprox
			<< endl;*/

		//cout << "\n  t = " << t << " <-> y: " << y.transpose() << endl;

		iter++;		
	} // while

	auto finish = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(finish - start);

	cout << "\nTime is " << 1e-6 * duration.count() << " sec.\n";

	cout << "\nNumber of steps: " << iter << endl;
	cout << "\n(t1 - t0) / hInitial = " << (t1 - t0) / hInitial << endl;
	cout << "\nMax step: " << hMax << endl;
	cout << "\nMax error approximated: " << errApprox << endl;

	cout << "\nFinal t: " << t << endl;
}
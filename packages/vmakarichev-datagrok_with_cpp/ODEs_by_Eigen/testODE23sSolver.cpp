// testODE23sSolver.cpp

// Implementation of the function testODE23sSolver()

/* Here, we test stiff equaltions solving method decribed in the papers
   https://doi.org/10.1137/S1064827594276424 and https://doi.org/10.1016/S0898-1221(00)00175-9
*/

#include<iostream>
#include<iomanip>
#include<ctime>
#include<fstream>
using namespace std;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "tests.h"

// STIFF EXAMPLE: Steve Chapra and Raymond P. Canale. Numerical Methods for Engineers, 2021 (page 770).

VectorXd F(double t, VectorXd & y)
{
	VectorXd res(y.size());

	res(0) = -5.0 * y(0) + 3.0 * y(1);
	res(1) = 100.0 * y(0) - 301.0 * y(1);

	return res;
}

MatrixXd J(double t, VectorXd & y)
{
	MatrixXd res(y.size(), y.size());

	res << -5.0, 3.0,
		100.0, -301.0;

	return res;
}

VectorXd T(double t, VectorXd & y)
{
	VectorXd res(y.size());

	res << 0.0, 
		0.0;

	return y;
}

VectorXd exactSolution(double t, unsigned dim)
{
	VectorXd res(dim);

	res(0) = 52.96 * exp(-3.9899 * t) - 0.67 * exp(-302.0101 * t);
	res(1) = 17.83 * exp(-3.9899 * t) + 65.99 * exp(-302.0101 * t);

	return res;
}


/*VectorXd F(double t, VectorXd & y)
{
	VectorXd res(y.size());

	res(0) = 2.0 * t;
	res(1) = 3.0 * t * t;

	return res;
}

MatrixXd J(double t, VectorXd & y)
{
	MatrixXd res(y.size(), y.size());

	res << 0.0, 0.0,
		0.0, 0.0;

	return res;
}

VectorXd T(double t, VectorXd & y)
{
	VectorXd res(y.size());

	res << 2.0,
		(6.0 * t);

	return y;
}

VectorXd exactSolution(double t, unsigned dim)
{
	VectorXd res(dim);

	res(0) = t * t;
	res(1) = t * t * t;

	return res;
} */

/*
VectorXd F(double t, VectorXd & y)
{
	VectorXd res(y.size());

	res(0) = 2.0 * exp(2.0 * t);
	res(1) = -3.0 * exp(-3.0 * t);

	return res;
}

MatrixXd J(double t, VectorXd & y)
{
	MatrixXd res(y.size(), y.size());

	res << 0.0, 0.0,
		0.0, 0.0;

	return res;
}

VectorXd T(double t, VectorXd & y)
{
	VectorXd res(y.size());

	res << (4.0 * exp(2.0 * t)),
		(9.0 * exp(-3.0 * t));

	return y;
}

VectorXd exactSolution(double t, unsigned dim)
{
	VectorXd res(dim);

	res(0) = exp(2.0 * t);
	res(1) = exp(-3.0 * t);

	return res;
}*/



// test ode23s solver
void testODE23sSolver()
{
	cout << "Analysis of solver for stiff problem: multi-dimensional case.\n";

	unsigned dim = 2;

	double h = 0.01;

	double t0 = 0.0;
	double t1 = 100.0;

	unsigned N = static_cast<unsigned>((t1 - t0) / h) + 1;

	cout << "\nThe initial problem is solved on the segment [t0, t1] with the step h,\nwhere\n"
		<< "   t0 = " << t0
		<< "\n   t1 = " << t1
		<< "\n   h = " << h
		<< "\n\n   number of points = " << N
		<< endl;

	double d = 1.0 - sqrt(2.0) / 2.0; // = 1 / (2 + sqrt(2))

	double t = t0;
		
	VectorXd y = exactSolution(t, dim);

	auto I = MatrixXd::Identity(dim, dim);

	MatrixXd W(dim, dim);

	MatrixXd invW(dim, dim);

	VectorXd f0(dim);

	VectorXd k1(dim);

	VectorXd f1(dim);

	VectorXd k2(dim);

	VectorXd yDer(dim);

	double mad = 0.0;

	double maxMad = 0.0;

	double maxRelative = 0.0;
	
	/*cout << "\ny: " << y.transpose() << endl;

	cout << "\nJ:\n" << J(t, y) << endl;

	cout << "\nT:\n" << T(t, y) << endl;
	
	cout << "\nd = " << d << endl; 

	cout << "\nI:\n" << I << endl; */

	ofstream file("ode32s.csv", ios::out);

	auto start = time(0);

	for (unsigned i = 0; i < N; i++)
	{
		f0 = h * F(t, y);

		W = I - h * d * J(t, y);

		invW = W.inverse();

		k1 = invW * (f0 + h * d * T(t, y));

		yDer = y + 0.5 * h * k1;

		f1 = F(t + 0.5 * h, yDer);

		k2 = invW * (f1 - k1) + k1;

		y += k2 * h;

		t += h;

		yDer = exactSolution(t, dim);

		mad = fabs(y(0) - yDer(0));

		for (unsigned i = 1; i < dim; i++)
			mad = fmax(mad, fabs(y(i) - yDer(i)));

		maxMad = fmax(mad, maxMad);

		file << mad << '\n';
		/*
		//cout << "\nmad = " << mad << endl;

		mad = fabs((y(0) - yDer(0)) / yDer(0));

		for (unsigned i = 1; i < dim; i++)
			mad = fmax(mad, fabs((y(i) - yDer(i)) / yDer(i)));

		maxRelative = fmax(mad, maxRelative);

		//cout << "relative = " << mad << endl;

		/*cout << "\niter no. " << i << endl;

		cout << "\nf0:\n" << f0 << endl;
		
		cout << "\nW:\n" << W << endl;		

		cout << "\ninverse W:\n" << invW << endl;

		cout << "\nW * invW:\n" << W * invW << endl;

		cout << "\nk1:\n" << k1 << endl;

		cout << "\nf1:\n" << f1 << endl;

		cout << "\nk2:\n" << k2 << endl;

		cout << "\nh * k2:\n" << h * k2 << endl;

		cout << "\ny: " << y.transpose() << endl;

		cout << "\nexact solution: " << yDer.transpose() << endl;*/

		//cout << "\n===============================================================\n";
	}

	auto finish = time(0);

	cout << "\nComputation time: " << finish - start << " sec.\n";

	cout << "\nMAXIMUM ABSOLUTE DEVIATION: " << maxMad << endl;
	cout << "\nMAXIMUM RELATIVE ERROR: " << maxRelative << endl;
	
}
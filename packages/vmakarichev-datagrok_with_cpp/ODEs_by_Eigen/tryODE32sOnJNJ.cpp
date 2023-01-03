// tryODE32sOnJNJ.cpp

// Implementation of the function tryODE32sOnJNJ()

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

namespace test {

	const int DIM = 2;

	VectorXd F(double t, VectorXd & y)
	{
		VectorXd res(y.size());

		res(0) = -5.0 * y(0) + 3.0 * y(1);
		res(1) = 100.0 * y(0) - 301.0 * y(1);

		//res(0) = y(0) * y(0) - y(1) * y(1) * y(1);
		//res(1) = 2.0 * y(0) + 13.0 * y(1);

		return res;
	}

	MatrixXd Japprox(double t, VectorXd & y, double eps)
	{
		MatrixXd res(y.size(), y.size());

		VectorXd val = F(t, y);		

		VectorXd yDer = y;

		for (int i = 0; i < y.size(); i++)
		{
			yDer(i) += eps;			

			res.col(i) = (F(t, yDer) - val) / eps;

			yDer(i) -= eps;
		}

		return res;
	}

	MatrixXd J(double t, VectorXd & y)
	{
		MatrixXd res(y.size(), y.size());

		res << -5.0, 3.0,
			100.0, -301.0;

		//res << (2.0 * y(0)), (-3.0 * y(1) * y(1)),
		//	2.0, 13.0;

		return res;
	}

	VectorXd Tapprox(double t, VectorXd & y, double eps)
	{
		return (F(t + eps, y) - F(t, y)) / eps;

		/*VectorXd res(y.size());

		res << 0.0,
			0.0;

		return y;*/
	}

	VectorXd exactSolution(double t, unsigned dim)
	{
		VectorXd res(dim);

		res(0) = 52.96 * exp(-3.9899 * t) - 0.67 * exp(-302.0101 * t);
		res(1) = 17.83 * exp(-3.9899 * t) + 65.99 * exp(-302.0101 * t);

		return res;
	}
}; // namespace test

namespace bigger {

	const int DIM = 7;

	VectorXd F(double t, VectorXd & y)
	{
		VectorXd res(y.size());

		res(0) = y(0) + 1;
		res(1) = y(1) + 2;
		res(2) = y(2) + t;
		res(3) = y(3) + t * t;
		res(4) = y(4) + exp(t);
		res(5) = y(5) + exp(-t);
		res(6) = y(6);

		return res;
	}

	MatrixXd Japprox(double t, VectorXd & y, double eps)
	{
		MatrixXd res(y.size(), y.size());

		VectorXd val = F(t, y);

		VectorXd yDer = y;

		for (int i = 0; i < y.size(); i++)
		{
			yDer(i) += eps;

			res.col(i) = (F(t, yDer) - val) / eps;

			yDer(i) -= eps;
		}

		return res;
	}
		
	VectorXd Tapprox(double t, VectorXd & y, double eps)
	{		
		return (F(t + eps, y) - F(t, y)) / eps;
	}

	VectorXd exactSolution(double t, unsigned dim)
	{
		VectorXd res(dim);

		res(0) = 1.0; //2.1;
		res(1) = 1.0; //1.2;
		res(2) = 1.0; //1.2;
		res(3) = 1.0; //2.1;
		res(4) = 1.0; //1.2;
		res(5) = 1.0; //1.2;
		res(6) = 1.0; //1.2;

		return res;
	}
}; // namespace bigger

namespace jnj {

	const int DIM = 7;

	VectorXd F(double t, VectorXd & y)
	{
		VectorXd res(y.size());

		// JNJ formulas are removed!
		return res;
	}

	MatrixXd Japprox(double t, VectorXd & y, double eps)
	{
		MatrixXd res(y.size(), y.size());

		VectorXd val = F(t, y);

		VectorXd yDer = y;

		for (int i = 0; i < y.size(); i++)
		{
			yDer(i) += eps;

			res.col(i) = (F(t, yDer) - val) / eps;

			yDer(i) -= eps;
		}

		return res;
	}

	VectorXd Tapprox(double t, VectorXd & y, double eps)
	{
		return (F(t + eps, y) - F(t, y)) / eps;
	}

	VectorXd exactSolution(double t, unsigned dim)
	{
		VectorXd res(dim);	
// JNJ formulas are removed!

		return res;
	}
}; // namespace jnj


// try ode32s solver on the jnj system
void tryODE32sOnJNJ()
{
	using namespace jnj;

	cout << "Analysis of solver for the JNJ problem.\n";

	unsigned dim = DIM;

	double h = 1e-3;
		//0.00000001;
	double eps = h / 10.0;

	double t0 = 0.0;
	double t1 = 1000.0;

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

	/*cout << "\nJ excact:\n" << J(t, y) << endl;

	cout << "\nJ approx:\n" << Japprox(t, y, h / 10.0) << endl;*/

	auto I = MatrixXd::Identity(dim, dim);

	MatrixXd W(dim, dim);

	MatrixXd invW(dim, dim);

	VectorXd f0(dim);

	VectorXd k1(dim);

	VectorXd f1(dim);

	VectorXd k2(dim);

	VectorXd yDer(dim);		

	ofstream file("ode32s_jnj.csv", ios::out);

	file << "t";

	for (unsigned i = 0; i < dim; i++)
		file << ",y_" << i << "(t)";

	file << endl;

	file << t;

	for (unsigned i = 0; i < dim; i++)
		file << "," << y(i);

	file << endl;

	auto start = time(0);

	double difference = 0.0;

	for (unsigned i = 0; i < N; i++)
	{
		f0 = h * F(t, y);

		//W = I - h * d * J(t, y); // usage of exact Jacobian

		W = I - h * d * Japprox(t, y, eps); // usage of approximated Jacobian

		invW = W.inverse();

		k1 = invW * (f0 + h * d * Tapprox(t, y, eps));

		yDer = y + 0.5 * h * k1;

		f1 = F(t + 0.5 * h, yDer);

		k2 = invW * (f1 - k1) + k1;

		yDer = y;

		y += k2 * h;

		difference = fmax(difference, ((y - yDer) / h - F(t, yDer)).cwiseAbs().maxCoeff());

		//cout << "\n  i = " << i << "   difference: " << ((y - yDer) / h - F(t,yDer) ).transpose() << endl;

		t += h;

		file << t;

		for (unsigned j = 0; j < dim; j++)
			file << "," << y(j);

		file << endl;
	}

	auto finish = time(0);

	cout << "\nComputation time: " << finish - start << " sec.\n";	

	cout << "\nMax difference: " << difference << endl;

}
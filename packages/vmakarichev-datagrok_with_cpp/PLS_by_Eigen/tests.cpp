#include<iostream>
#include<cstdlib>
#include<ctime>
using namespace std;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

//typedef float Float;

#include "tests.h"

#include "PLS.h"
using namespace pls;

void tryPLS1_ver1()
{
	Float inputData[] = {1, 22, 3, 4, 25, 16, 8, 10, 121, 14, 16, 329, 22, 25, 208};

	const int ROWS_COUNT = 5;
	const int COLUMNS_COUNT = 3;
	const int L = COLUMNS_COUNT;// 1;

	Float outputData[ROWS_COUNT] =
		//{ -1, 0, 1, 0, -1 };
		// { 1, -2, 0, -4, 5};
	{ 0 };

	Map < Matrix<Float, Dynamic, Dynamic, ColMajor>> D(inputData, ROWS_COUNT, COLUMNS_COUNT);

	Map < Vector<Float, Dynamic>> y(outputData, ROWS_COUNT);

	for (int i = 0; i < COLUMNS_COUNT; i++)
		y += (i + 1) * D.col(i);

	cout << "y:\n" << y << endl;

	cout << "\nD:\n" << D << endl;

	Vector<Float, Dynamic> mu = D.colwise().mean();

	cout << "\nmu:\n" << mu << endl;

	Matrix<Float, Dynamic, Dynamic, ColMajor> X = D.rowwise() - mu.transpose();


	
	Matrix<Float, Dynamic, Dynamic, ColMajor> W(COLUMNS_COUNT, L);

	Matrix<Float, Dynamic, Dynamic, ColMajor> P(COLUMNS_COUNT, L);

	Vector<Float, Dynamic> q(L);

	for (int k = 0; k < L; k++)
	{
		Vector<Float, Dynamic> w = (X.transpose() * y).normalized();

		W.col(k) = w;

		Vector<Float, Dynamic> t = (X * w).normalized();

		Vector<Float, Dynamic> p = X.transpose() * t;

		P.col(k) = p;

		X = X - t * (p.transpose());

		q(k) = t.transpose() * y;
	}

	cout << "\nW:\n" << W << endl;

	cout << "\nP:\n" << P << endl;

	Vector<Float, Dynamic> b = W * (P.transpose() * W).inverse() * q;

	cout << "\nb:\n" << b << endl;

	cout << "\nq:\n" << q << endl;

	cout << "\nD:\n" << D << endl;

	cout << "\nD * b:\n" << D * b << endl;// (D.rowwise() - mu.transpose()) * b << endl;

	cout << "\ny:\n" << y << endl;
} // tryPLS1_ver1

void tryPLS1_ver2()
{
	Float inputData[] = { 1, 22, 3, 4, 25, 16, 8, 10, 121, 14, 16, 329, 22, 25, 208 };

	const int ROWS_COUNT = 5;
	const int COLUMNS_COUNT = 3;
	const int L = COLUMNS_COUNT;// 1;

	Float outputData[ROWS_COUNT] =
		//{ -1, 0, 1, 0, -1 };
		// { 1, -2, 0, -4, 5};
	{ 0 };

	Map < Matrix<Float, Dynamic, Dynamic, ColMajor>> D(inputData, ROWS_COUNT, COLUMNS_COUNT);

	Map < Vector<Float, Dynamic>> y(outputData, ROWS_COUNT);

	for (int i = 0; i < COLUMNS_COUNT; i++)
		y += (i + 1) * D.col(i);

	cout << "y:\n" << y << endl;

	cout << "\nD:\n" << D << endl;

	Vector<Float, Dynamic> mu = D.colwise().mean();

	cout << "\nmu:\n" << mu << endl;

	Matrix<Float, Dynamic, Dynamic, ColMajor> X = D.rowwise() - mu.transpose();
	
	Matrix<Float, Dynamic, Dynamic, ColMajor> W(COLUMNS_COUNT, L);

	Matrix<Float, Dynamic, Dynamic, ColMajor> P(COLUMNS_COUNT, L);

	Matrix<Float, Dynamic, Dynamic, ColMajor> T(ROWS_COUNT, L);

	Vector<Float, Dynamic> normTau(L);

	Vector<Float, Dynamic> q(L);

	Vector<Float, Dynamic> normV(L);


	Vector<Float, Dynamic> w = (X.transpose() * y);

	normV(0) = w.norm();

	w = w / normV(0);

	W.col(0) = w;

	Vector<Float, Dynamic> t = X * w;

	normTau(0) = t.norm();

	t = t / normTau(0);

	T.col(0) = t;

	Vector<Float, Dynamic> p = X.transpose() * t;

	P.col(0) = p;

	q(0) = t.transpose() * y;


	for (int a = 1; a < L; a++)
	{
		w = normV(a - 1) * (w - p / normTau(a - 1));

		normV(a) = w.norm();

		w = w / normV(a);

		W.col(a) = w;

		t = X * w;

		t = t - T.leftCols(a) * (T.leftCols(a).transpose() * t);

		normTau(a) = t.norm();

		t = t / normTau(a);

		T.col(a) = t;

		p = X.transpose() * t;

		P.col(a) = p;

		q(a) = t.transpose() * y;
	}
	

	cout << "\nW:\n" << W << endl;

	cout << "\nP:\n" << P << endl;

	Vector<Float, Dynamic> b = W * (P.transpose() * W).inverse() * q;

	cout << "\nb:\n" << b << endl;

	cout << "\nq:\n" << q << endl;

	cout << "\nD:\n" << D << endl;

	cout << "\nD * b:\n" << D * b << endl;// (D.rowwise() - mu.transpose()) * b << endl;

	cout << "\ny:\n" << y << endl;
} // tryPLS1_ver2

void testPLS1()
{
	Float inputData[] = { 1, 22, 3, 4, 25, 16, 8, 10, 121, 14, 16, 329, 22, 25, 208 };

	const int ROWS_COUNT = 5;
	const int COLUMNS_COUNT = 3;
	const int L = COLUMNS_COUNT;// 1;

	Float outputData[ROWS_COUNT]
		= { 0 };
		// = { -1, 0, 1, 0, -1 };
		// = { 1, -2, 0, -4, 5};

	Float prediction[ROWS_COUNT] = { 0 };

	cout << "X & Y:\n";
	for (int i = 0; i < ROWS_COUNT; i++) 
	{
		for (int j = 0; j < COLUMNS_COUNT; j++)
		{
			outputData[i] += (j + 1) * inputData[i + j * ROWS_COUNT]; 
			cout << "  " << inputData[i + j * ROWS_COUNT];
		}
		cout << "     " << outputData[i] << endl;
	}

	Float coefs[ROWS_COUNT] = { 0 };

	cout << "PLS result code: " 
		<< partialLeastSquare(inputData, ROWS_COUNT, COLUMNS_COUNT, outputData, L, prediction, coefs) << endl;

	cout << "\ncoefficients:\n";

	for (int i = 0; i < COLUMNS_COUNT; i++)
		cout << "  " << coefs[i];

	cout << "\n\nY  & its prediction:\n";
	for (int i = 0; i < ROWS_COUNT; i++)
		cout << "  " << outputData[i] << "   " << prediction[i] << endl;

	cout << "\nMaximum absolute deviation between Y & its prediction: "
		<< mad(outputData, prediction, COLUMNS_COUNT) << endl;
	
} // testPLS1

void plsPerformance()
{		
	const int SEED = time(0);//342311;
	const int ROWS_COUNT = 1000000;
	const int COLUMNS_COUNT = 150;
	const int COMPONENTS_COUNT = 3;// COLUMNS_COUNT;
	const int ITERATIONS_NUM = 20;

	srand(SEED);

	Float * inputData = new Float[ROWS_COUNT * COLUMNS_COUNT];

	for (int i = 0; i < ROWS_COUNT * COLUMNS_COUNT; i++)
		inputData[i] = static_cast<Float>( rand() % 100 );

	Float * linCombCoef = new Float[COLUMNS_COUNT];

	for (int i = 0; i < COLUMNS_COUNT; i++)
		linCombCoef[i] = static_cast<Float>(-10 + rand() % 21);//(i + 1);

	/*cout << "\ninit coefficients:\n";
	for (int i = 0; i < COLUMNS_COUNT; i++)
		cout << "  " << linCombCoef[i];*/

	Float * outputData = new Float[ROWS_COUNT];

	Float* prediction = new Float[ROWS_COUNT];
		
	//cout << "\nX & Y:\n";
	for (int i = 0; i < ROWS_COUNT; i++)
	{
		outputData[i] = 0;
		prediction[i] = 0;

		for (int j = 0; j < COLUMNS_COUNT; j++)
		{
			outputData[i] += linCombCoef[j] * inputData[i + j * ROWS_COUNT];
			//cout << "  " << inputData[i + j * ROWS_COUNT];
		}
		//cout << "     " << outputData[i] << endl;
	}

	Float * coefs = new Float[COLUMNS_COUNT];

	Float sum = 0;
	time_t start = 0;
	time_t finish = 0;

	for (int i = 1; i <= ITERATIONS_NUM; i++)
	{
		start = time(0);

		cout << "Launch " << i 
			<<"\n    result code: "
			<< partialLeastSquare(inputData, ROWS_COUNT, COLUMNS_COUNT, outputData, COMPONENTS_COUNT, prediction, coefs) 
			<< "\n";

		finish = time(0);

		cout << "    time: " << finish - start << " sec.\n\n";

		sum += finish - start;
	}

	cout << "\nAverage time: " << sum / ITERATIONS_NUM << endl;

	/*cout << "\ninit coefficients:\n";
	for (int i = 0; i < COLUMNS_COUNT; i++)
		cout << "  " << linCombCoef[i];
	
	cout << "\ncomputed coefficients:\n";
	for (int i = 0; i < COLUMNS_COUNT; i++)
		cout << "  " << coefs[i];*/

	cout << "\nMaximum absolute deviation between initial and computed coefficients: "
		<< mad(linCombCoef, coefs, COLUMNS_COUNT) << endl;

	/*cout << "\n\nY  & its prediction:\n";
	for (int i = 0; i < ROWS_COUNT; i++)
		cout << "  " << outputData[i] << "   " << prediction[i] << endl;*/

	cout << "\nMaximum absolute deviation between Y & its prediction: "
		<< mad(outputData, prediction, ROWS_COUNT) << endl;

	delete[] inputData;
	delete[] outputData;
	delete[] prediction;
	delete[] linCombCoef;
	delete[] coefs;
} // plsPerformance

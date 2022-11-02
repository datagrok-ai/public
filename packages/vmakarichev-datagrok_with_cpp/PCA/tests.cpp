// tests.cpp
// Implementations of test-functions

#include<iostream>
#include<cstdlib>
#include<ctime>

using namespace std;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "tests.h"
#include "PCA.h"
#include "constants.h"
using namespace pca;

// Manipulating matrices of different types: PROBLEM when different types!
void differentTypeValuesManipulation()
{
	// put 'int' or 'float' instead of 'double', and an error occurs!
	Matrix<double, 2, 2> A;
	A << 1, 2,
		3, 4;

	cout << "A:\n" << A << endl;

	Matrix<double, 2, 2> B;
	B << 1.1, 2.2,
		3.3, 4.4;

	cout << "\nB:\n" << B << endl;

	cout << "\nA + B:\n" << A + B << endl;

	cout << "\nA * B:\n" << A * B << endl;
} // differentTypeValuesManipulation

// Test of PCA that is performed using correlation matrix
void pcaUsingCorrelationMatrixTest()
{
	const int width = 6;

	// float data
	const int heightOfFloats = 3;
	float a[width] = {103.1f, 202.2f, 305.3f, 410.4f, 505.5f, 606.6f};
	float b[width] = {-106.1f, -216.2f, -390.3f, -441.4f, -552.5f, -606.f};
	float c[width] = { 103.1f, -209.2f, 380.3f, -400.4f, 570.5f, 600.7f};

	// int data
	const int heightOfInts = 2;
	int x[width] = { 64, 315, 205, 114, 196, 178 };
	int y[width] = { 152, 61, 171, 212, 121, 290 };	

	const int height = heightOfFloats + heightOfInts;

	// test data
	void ** data = new void * [height];
	data[0] = a;
	data[1] = b;
	data[2] = c;
	data[3] = x;
	data[4] = y;

	cout << "Test data:\n";

	for (int i = 0; i < heightOfFloats; i++)
	{
		float * ptr = static_cast<float *>(data[i]);
		for (int j = 0; j < width; j++)
			cout << ptr[j] << "  ";
		cout << endl;
	} 

	for (int i = heightOfFloats; i < height; i++)
	{
		int * ptr = static_cast<int *>(data[i]);
		for (int j = 0; j < width; j++)
			cout << ptr[j] << "  ";
		cout << endl;
	}

	
	const int numOfPrincipalComponents = 3;

	float ** principalComponents = new float* [numOfPrincipalComponents];
	for (int i = 0; i < numOfPrincipalComponents; i++)
		principalComponents[i] = new float [width];

	// compute principal components
	pcaUsingCorrelationMatrix(data, heightOfFloats, heightOfInts, width, numOfPrincipalComponents, principalComponents);

	cout << "\nPrincipal components (data is given by void **):\n";
	for (int i = 0; i < numOfPrincipalComponents; i++)
	{
		for (int j = 0; j < width; j++)
			cout << "  " << principalComponents[i][j];
		cout << endl;
	}
	
	// Now, we compute principal components using the second approach: data is given by float *

	float matrix[] = { 103.1f, 202.2f, 305.3f, 410.4f, 505.5f, 606.6f,
		-106.1f, -216.2f, -390.3f, -441.4f, -552.5f, -606.f,
		103.1f, -209.2f, 380.3f, -400.4f, 570.5f, 600.7f,
		64, 315, 205, 114, 196, 178,
		152, 61, 171, 212, 121, 290 };

	float princComp[numOfPrincipalComponents * width];

	pcaUsingCorrelationMatrix(matrix, height, width, numOfPrincipalComponents, princComp);

	cout << "\nPrincipal components (data is given by float *):\n";

	for (int i = 0; i < numOfPrincipalComponents; i++)
	{
		for (int j = 0; j < width; j++)
			cout << "  " << princComp[i * width + j];
		cout << endl;
	}

	/*
		
	Vector<float, Dynamic> mu(heightOfFloats + heightOfInts);

	mu = D.rowwise().mean();

	Matrix<float, Dynamic, Dynamic> C(heightOfFloats + heightOfInts, heightOfFloats + heightOfInts);

	C = D * D.transpose() / width - mu * mu.transpose();

	SelfAdjointEigenSolver<Matrix<float, Dynamic, Dynamic>> eigensolver(C);
	
	Matrix<float, Dynamic, Dynamic, ColMajor> P = eigensolver.eigenvectors().rowwise().reverse();

	Matrix<float, Dynamic, Dynamic> E = P(all, seq(0, numOfPrincipalComponents - 1));
	
	//Matrix<float, Dynamic, Dynamic> Princ = (D.colwise() - mu).transpose() * E;

	Matrix<float, Dynamic, Dynamic> Princ = E.transpose() * (D.colwise() - mu) ;

	cout << "\nPrincipal components:\n" << Princ  << endl;

	//cout << "\nApproximation:\n" << (E * Princ.transpose()).colwise() + mu << endl;*/
		
	delete [] data;
	
	for (int i = 0; i < numOfPrincipalComponents; i++)
		delete[] principalComponents[i];

	delete[] principalComponents;
} // pcaUsingCorrelationMatrixTest

// Performance comparison of PCAs implemented using correlation matrix
void comparePerformanceOfPCAwithCorMatr()
{	
	// set random seed
	srand(SEED);
	
	// Set test data
	void ** data = new void *[HEIGHT];

	for (int i = 0; i < HEIGHT_OF_FLOATS; i++)
	{
		float * fPtr = new float [WIDTH];
		for (int j = 0; j < WIDTH; j++)
			fPtr[j] = static_cast<float>(rand() % 200) / (1 + rand() % 10);		
		data[i] = fPtr;
	}
	for (int i = HEIGHT_OF_FLOATS; i < HEIGHT; i++)
	{
		int * iPtr = new int [WIDTH];
		for (int j = 0; j < WIDTH; j++)
			iPtr[j] = rand() % 20;
		data[i] = iPtr;
	}	
	
	/*cout << "Test data:\n";
	for (int i = 0; i < heightOfFloats; i++)
	{
		float * ptr = static_cast<float *>(data[i]);
		for (int j = 0; j < width; j++)
			cout << ptr[j] << "  ";
		cout << endl;
	}

	for (int i = heightOfFloats; i < height; i++)
	{
		int * ptr = static_cast<int *>(data[i]);
		for (int j = 0; j < width; j++)
			cout << ptr[j] << "  ";
		cout << endl;
	}*/
	
	// Principal components: float **
	float ** principalComponents = new float *[NUM_OF_PRINCIPAL_COMPONENTS];
	for (int i = 0; i < NUM_OF_PRINCIPAL_COMPONENTS; i++)
		principalComponents[i] = new float[WIDTH];

	// compute principal components: data is given by void **
	time_t t0 = time(0);

	pcaUsingCorrelationMatrix(data, HEIGHT_OF_FLOATS, HEIGHT_OF_INTS, WIDTH, 
		NUM_OF_PRINCIPAL_COMPONENTS, principalComponents);

	time_t t1 = time(0);

	cout << "Time of PCA (data is given by void **):  " << t1 - t0 << " sec.\n";

	// Now, we compute principal components using the second approach: data is given by float *

	float * matrix = new float[HEIGHT * WIDTH];

	// set data: float part
	for (int i = 0; i < HEIGHT_OF_FLOATS; i++)
	{
		float * fPtr = static_cast<float *>(data[i]);
		for (int j = 0; j < WIDTH; j++)
			matrix[i * WIDTH + j] = fPtr[j];
	}
	// set data: int part
	for (int i = HEIGHT_OF_FLOATS; i < HEIGHT; i++)
	{
		int * iPtr = static_cast<int *>(data[i]);
		for (int j = 0; j < WIDTH; j++)
			matrix[i * WIDTH + j] = static_cast<float>(iPtr[j]);
	}

	// principal components
	float * princComp = new float [NUM_OF_PRINCIPAL_COMPONENTS * WIDTH];

	// compute principal components: data is given by float *
	t0 = time(0);
	
	pcaUsingCorrelationMatrix(matrix, HEIGHT, WIDTH, NUM_OF_PRINCIPAL_COMPONENTS, princComp);	

	t1 = time(0);

	cout << "\nTime of PCA (data is given by float *):  " << t1 - t0 << " sec.\n";
		
	// Compute maximum absolute deviation (MAD) between the results obtained using different approaches
	float mad = 0.0;

	for (int i = 0; i < NUM_OF_PRINCIPAL_COMPONENTS; i++)
		for (int j = 0; j < WIDTH; j++)
			mad = max(mad, fabsf(principalComponents[i][j] - princComp[i * WIDTH + j]));

	cout << "\nMaximum absolute deviation between obtained results:  " << mad << endl;

	// Clear allocated memory

	for (int i = 0; i < HEIGHT; i++)
		delete[] data[i];

	delete[] data;

	for (int i = 0; i < NUM_OF_PRINCIPAL_COMPONENTS; i++)
		delete[] principalComponents[i];

	delete[] principalComponents;

	delete[] matrix;

	delete[] princComp;
} // comparePerformanceOfPCAwithCorMatr

// Performance investiagtion of PCA using correlation matrix: data is given by void **
void investigatePCAwithCorMatrVoidDataRepresentation(int numOfLaunches)
{
	// set random seed
	srand(SEED);

	// Set test data
	void ** data = new void *[HEIGHT];

	// set int part
	for (int i = 0; i < HEIGHT_OF_FLOATS; i++)
	{
		float * fPtr = new float[WIDTH];
		for (int j = 0; j < WIDTH; j++)
			fPtr[j] = static_cast<float>(rand() % 200) / (1 + rand() % 10);
		data[i] = fPtr;
	}

	// set float part
	for (int i = HEIGHT_OF_FLOATS; i < HEIGHT; i++)
	{
		int * iPtr = new int[WIDTH];
		for (int j = 0; j < WIDTH; j++)
			iPtr[j] = rand() % 20;
		data[i] = iPtr;
	}

	// Principal components: float **
	float ** principalComponents = new float *[NUM_OF_PRINCIPAL_COMPONENTS];
	for (int i = 0; i < NUM_OF_PRINCIPAL_COMPONENTS; i++)
		principalComponents[i] = new float[WIDTH];
	
	// vector to store time
	Vector<float, Dynamic> numOfSecondsPerLaunch(numOfLaunches);

	cout << "Time for PCA (data is given by void **).\n";

	// Launch computing process
	for (int i = 1; i <= numOfLaunches; i++)
	{
		time_t t0 = time(0);

		pcaUsingCorrelationMatrix(data, HEIGHT_OF_FLOATS, HEIGHT_OF_INTS, WIDTH,
			NUM_OF_PRINCIPAL_COMPONENTS, principalComponents);

		time_t t1 = time(0);

		cout << "  Launch " << i << ":  " << t1 - t0 << " sec.\n";

		numOfSecondsPerLaunch(i - 1) = static_cast<float>(t1 - t0);
	}

	cout << "\nAverage time per launch: " << numOfSecondsPerLaunch.mean() << " sec.";

	// Clear allocated memory

	for (int i = 0; i < HEIGHT; i++)
		delete[] data[i];

	delete[] data;

	for (int i = 0; i < NUM_OF_PRINCIPAL_COMPONENTS; i++)
		delete[] principalComponents[i];

	delete[] principalComponents;
} // investigatePCAwithCorMatrVoidDataRepresentation

// Performance investiagtion of PCA using correlation matrix: data is given by float *
void investigatePCAwithCorMatrFloatDataRepresentation(int numOfLaunches)
{
	// set random seed
	srand(SEED);

	// Set test data

	float * matrix = new float[HEIGHT * WIDTH];

	for (int i = 0; i < HEIGHT_OF_FLOATS; i++)		
		for (int j = 0; j < WIDTH; j++)
			matrix[i * WIDTH + j] = static_cast<float>(rand() % 200) / (1 + rand() % 10);

	for (int i = HEIGHT_OF_FLOATS; i < HEIGHT; i++)
		for (int j = 0; j < WIDTH; j++)
			matrix[i * WIDTH + j] = static_cast<float>(rand() % 20);
		
	// principal components
	float * princComp = new float[NUM_OF_PRINCIPAL_COMPONENTS * WIDTH];

	// vector to store time
	Vector<float, Dynamic> numOfSecondsPerLaunch(numOfLaunches);

	cout << "Time for PCA (data is given by float *).\n";

	// Launch computing process
	for (int i = 1; i <= numOfLaunches; i++)
	{
		time_t t0 = time(0);

		pcaUsingCorrelationMatrix(matrix, HEIGHT, WIDTH, NUM_OF_PRINCIPAL_COMPONENTS, princComp);

		time_t t1 = time(0);

		cout << "  Launch " << i << ":  " << t1 - t0 << " sec.\n";

		numOfSecondsPerLaunch(i - 1) = static_cast<float>(t1 - t0);
	}

	cout << "\nAverage time per launch: " << numOfSecondsPerLaunch.mean() << " sec.";
	
	delete[] matrix;

	delete[] princComp;
} // comparePerformanceOfPCAwithCorMatr
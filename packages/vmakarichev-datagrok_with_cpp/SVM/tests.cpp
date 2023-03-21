// tests.cpp

// Implementations of funtions for testing methods.

#include<iostream>
#include<cstring>
#include<fstream>
using namespace std;

#include "svm.h"
using namespace svm;

#include "generators.h"
using namespace gener;

#include "tests.h"

void testCreateDataset()
{
	const int ROW_COUNT = 5;
	const int COL_COUNT = 3;
	const int SIZE = ROW_COUNT * COL_COUNT;

	float data[SIZE] = { 0 };
	
	int index = 0;

	for (int j = 0; j < COL_COUNT; j++)
	{
		for (int i = 0; i < ROW_COUNT; i++)
			data[index++] = (float)j;
	}

	cout << "Data:\n";
	for(int i = 0; i < SIZE; i++)
		cout << "  " << data[i];

	float dataset[SIZE] = { 0 };

	cout << "\nDataset before:\n";
	for (int i = 0; i < SIZE; i++)
		cout << "  " << dataset[i];
	
	createDataset(data, ROW_COUNT, COL_COUNT, dataset);

	cout << "\nDataset after:\n";
	for (int i = 0; i < SIZE; i++)
		cout << "  " << dataset[i];

} // testCreateDataset

void testTrainModelSimpleLinear()
{
	const int ROW_COUNT = 6;
	const int COL_COUNT = 2;
	const int SIZE = ROW_COUNT * COL_COUNT;

	float data[SIZE] = { 1.0f, 2.0f, 0.5, -1.0f, -2.0f, -0.5, 2.0f, 1.0f, 0.5f, -2.0f, -1.0f, -0.5f };

	cout << "Data:\n";
	for (int i = 0; i < SIZE; i++)
		cout << "  " << data[i];

	float xTrain[SIZE] = { 0 };

	createDataset(data, ROW_COUNT, COL_COUNT, xTrain);

	cout << "\nX train:\n";
	for (int i = 0; i < SIZE; i++)
		cout << "  " << xTrain[i];

	float yTrain[ROW_COUNT] = {1.0f, 1.0f, 1.0f, -1.0f, -1.0f, -1.0f };

	float gamma = 1.0f;
	int kernel = LINEAR;
	float kernelParams[MAX_NUM_OF_KERNEL_PARAM] = {1.0f, 1.0f};

	float modelParams[ROW_COUNT + 1] = { 0.0f };

	cout << "\n\nTrain model result: "
		<< trainLSSVM(gamma, kernel, kernelParams, xTrain, yTrain, ROW_COUNT, COL_COUNT, modelParams)
		<< endl;

	cout << "\nParameters are:\n";
	for (int i = 0; i <= ROW_COUNT; i++)
		cout << "  " << modelParams[i];
} // testTrainModelSimpleLinear

void testGeneratorLinearSeparable()
{
	const int featuresCount = 2;
	const int samplesCount = 100000;

	float w[featuresCount] = {-1, 3};
	float b = -1;

	float minVal[featuresCount] = { -5, -2 };
	float maxVal[featuresCount] = { 5, 7 };

	string fileName = "svm_dataset_-1_3_-1_100000_2.csv";

	float * data = new float [featuresCount * samplesCount];
	float * labels = new float [samplesCount];

	cout << "Generate data: "
		<< generateLinearSeparable(w, b, featuresCount, samplesCount,
			minVal, maxVal, data, labels) << endl;

	/*cout << "\ndata:\n";
	for (int i = 0; i < featuresCount * samplesCount; i++)
		cout << "  " << data[i];
	
	cout << "\nlabels:\n";
	for (int i = 0; i < samplesCount; i++)
		cout << "  " << labels[i];*/

	ofstream file(fileName, ios::out);
	if (!file)
	{
		cout << "FAIL TO CREATE FILE FOR RESULTS!\n";
		return;
	}

	for (int j = 1; j <= featuresCount; j++)
		file << 'x' << j << ',';
	file << "y\n";

	for (int i = 0; i < samplesCount; i++)
	{
		for (int j = 0; j < featuresCount; j++)
			file << data[i + j * samplesCount] << ',';

		file << labels[i] << endl;
	}

	delete[] data;
	delete[] labels;

} // testGeneratorLinearSeparable

void testGeneratorLinearNonSeparable()
{
	const int featuresCount = 2;
	const int samplesCount = 10000;

	float w[featuresCount] = { -3, -2 };
	float b = 3;

	float minVal[featuresCount] = { -5, -2 };
	float maxVal[featuresCount] = { 5, 7 };

	float violatorsPercentage = 10;
	float actualViolatorsPercentage = 0;

	string fileName = "svm_dataset_-3_-2_3_10000_2_Non.csv";

	float* data = new float[featuresCount * samplesCount];
	float* labels = new float[samplesCount];

	cout << "Generate data: "
		<< generateLinearNonSeparable(w, b, featuresCount, samplesCount,
			minVal, maxVal, data, labels, violatorsPercentage, actualViolatorsPercentage) << endl;

	cout << "\nActual violators: " << actualViolatorsPercentage << endl;

	/*cout << "\ndata:\n";
	for (int i = 0; i < featuresCount * samplesCount; i++)
		cout << "  " << data[i];

	cout << "\nlabels:\n";
	for (int i = 0; i < samplesCount; i++)
		cout << "  " << labels[i];*/

	//cout << "\nChanging labels: " << changeLabels(labels, samplesCount, changeProbability) << endl;

	/*cout << "\nlabels after changing:\n";
	for (int i = 0; i < samplesCount; i++)
		cout << "  " << labels[i];*/

	ofstream file(fileName, ios::out);
	if (!file)
	{
		cout << "FAIL TO CREATE FILE FOR RESULTS!\n";
		return;
	}

	for (int j = 1; j <= featuresCount; j++)
		file << 'x' << j << ',';
	file << "y\n";

	for (int i = 0; i < samplesCount; i++)
	{
		for (int j = 0; j < featuresCount; j++)
			file << data[i + j * samplesCount] << ',';

		file << labels[i] << endl;
	}

	delete[] data;
	delete[] labels;

} // testGeneratorLinearNonSeparable

void testTrainModelComplexLinear()
{
	const int featuresCount = 2;
	const int samplesCount = 2000;

	float w[featuresCount] = { -1, 1 };
	float b = 0;

	/*float minVal[featuresCount] = { -5, -2 };
	float maxVal[featuresCount] = { 5, 7 };*/

	float minVal[featuresCount] = { -1, -1 };
	float maxVal[featuresCount] = { 1, 1 };

	float violatorsPercentage = 10;
	float actualViolatorsPercentage = 0;

	string fileName = "svm_test.csv";

	float* data = new float[featuresCount * samplesCount];
	float* labels = new float[samplesCount];
	float* xTrain = new float[featuresCount * samplesCount];
	float* yTrain = labels;

	cout << "Generate data: "
		<< generateLinearSeparable(w, b, featuresCount, samplesCount,
			minVal, maxVal, data, labels) << endl;

	/*cout << "Generate data: "
		<< generateLinearNonSeparable(w, b, featuresCount, samplesCount,
			minVal, maxVal, data, labels, violatorsPercentage, actualViolatorsPercentage) << endl;*/

	cout << "\nActual violators, %: " << actualViolatorsPercentage << endl;

	cout << "\nCreate xTrain: " << createDataset(data, samplesCount, featuresCount, xTrain) << endl;

	/*cout << "\nX train:\n";
	for (int i = 0; i < featuresCount * samplesCount; i++)
		cout << "  " << xTrain[i];

	cout << "\nyTrain:\n";
	for (int i = 0; i < samplesCount; i++)
		cout << "  " << yTrain[i];*/


	float gamma = 1.0f;
	int kernel = LINEAR;
	float kernelParams[MAX_NUM_OF_KERNEL_PARAM] = { 1.0f, 1.0f };

	float modelParams[samplesCount + 1] = { 0.0f };

	cout << "\n\nTrain model result: "
		<< trainLSSVM(gamma, kernel, kernelParams, xTrain, yTrain, samplesCount, featuresCount, modelParams)
		<< endl;

	/*cout << "\nParameters are:\n";
	for (int i = 0; i <= samplesCount; i++)
		cout << "  " << modelParams[i];*/

	cout << "\nbias: " << modelParams[samplesCount] << endl;


	ofstream file(fileName, ios::out);
	if (!file)
	{
		cout << "FAIL TO CREATE FILE FOR RESULTS!\n";
		return;
	}

	for (int j = 1; j <= featuresCount; j++)
		file << 'x' << j << ',';
	file << "y\n";

	for (int i = 0; i < samplesCount; i++)
	{
		for (int j = 0; j < featuresCount; j++)
			file << data[i + j * samplesCount] << ',';

		file << labels[i] << endl;
	}

	delete[] data;
	delete[] labels;
	delete[] xTrain;

} // testTrainModelComplexLinear
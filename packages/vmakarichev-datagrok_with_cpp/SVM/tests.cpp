// tests.cpp

// Implementations of funtions for testing methods.

#include<iostream>
using namespace std;

#include "svm.h"
using namespace svm;

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

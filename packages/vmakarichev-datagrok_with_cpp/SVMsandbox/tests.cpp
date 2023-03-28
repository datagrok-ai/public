#include<iostream>
#include<fstream>
using namespace std;

#include "svm.h"
using namespace svm;

#include "dataMining.h"
using dmt::getNormalizedDataset;

#include "svmApi.h"

#include "tests.h"

void testGenerateDataSet()
{
	int kernel = LINEAR;
	int kernelParamsCount = 2;
	float * kernelParams = new float [kernelParamsCount];
	int rowCount = 10000;
	int colCount = 2;
	int labelsLength = rowCount;
	float* dataset = new float[rowCount * colCount];
	float* labels = new float[labelsLength];
	float min = -153;
	float max = 791;
	float violatorsPercentage = 1;

	// RBF test
	kernel = RBF;
	float sigma = (max - min) / 4;
	kernelParams[0] = sigma;

	cout << "Generating dataset: "
		<< generateDataset(kernel, kernelParams, kernelParamsCount, 
			dataset, rowCount, colCount, labels, labelsLength, 
			min, max, violatorsPercentage)
		<< endl;

	string fileName = "dataset.csv";

	ofstream file(fileName, ios::out);
	if (!file)
	{
		cout << "FAIL TO CREATE FILE FOR RESULTS!\n";
		return;
	}

	for (int j = 1; j <= colCount; j++)
		file << 'x' << j << ',';
	file << "y\n";

	for (int i = 0; i < rowCount; i++)
	{
		for (int j = 0; j < colCount; j++)
			file << dataset[i + j * rowCount] << ',';

		file << labels[i] << endl;
	}

	delete[] kernelParams;
	delete[] dataset;
	delete[] labels;
}

void testNormalizeDataset()
{
	int kernel = LINEAR;
	int kernelParamsCount = 2;
	float* kernelParams = new float[kernelParamsCount];
	int rowCount = 10;
	int colCount = 2;
	int labelsLength = rowCount;
	float* dataset = new float[rowCount * colCount];
	float* normalizedDataset = new float[rowCount * colCount];
	float* means = new float[colCount];
	float* stdDevs = new float[colCount];
	float* labels = new float[labelsLength];
	float min = -10;
	float max = 10;
	float violatorsPercentage = 1;

	// RBF test
	kernel = RBF;
	float sigma = (max - min) / 4;
	kernelParams[0] = sigma;

	cout << "Generating dataset: "
		<< generateDataset(kernel, kernelParams, kernelParamsCount,
			dataset, rowCount, colCount, labels, labelsLength,
			min, max, violatorsPercentage)
		<< endl;
	
	cout << "Normalizing: "
		<< normalizeDataset(dataset, rowCount, colCount,
			normalizedDataset, colCount, rowCount,
			means, colCount, stdDevs, colCount)
		<< endl;

	delete[] kernelParams;
	delete[] dataset;
	delete[] labels;
	delete[] normalizedDataset;
	delete[] means;
	delete[] stdDevs;
}

void testTrainLSSVM()
{
	float gamma = 1.0f;

	// Kernel
	int kernel = LINEAR;
	int kernelParamsCount = 2;
	float* kernelParams = new float[kernelParamsCount];

	// Data
	int rowCount = 600;
	int colCount = 2;
	int labelsLength = rowCount;
	float* dataset = new float[rowCount * colCount];
	float* normalizedDataset = new float[rowCount * colCount];
	float* means = new float[colCount];
	float* stdDevs = new float[colCount];
	float* labels = new float[labelsLength];
	float* predictedLabels = new float[labelsLength];
	float min = -10;
	float max = 10;
	float violatorsPercentage = 1;

	// Model
	int modelParamsCount = rowCount + 1;
	float* modelParams = new float[modelParamsCount];
	int precomputedWeightsCount = colCount + 1;
	float* precomputedWeights = new float[precomputedWeightsCount];
			
	cout << "Generating dataset: "
		<< generateDataset(kernel, kernelParams, kernelParamsCount,
			dataset, rowCount, colCount, labels, labelsLength,
			min, max, violatorsPercentage)
		<< endl;

	cout << "\nTraining model: "
		<< trainLSSVM(gamma, kernel, kernelParams, kernelParamsCount,
			modelParamsCount, precomputedWeightsCount,
			dataset, rowCount, colCount, labels, labelsLength,
			normalizedDataset, colCount, rowCount,
			means, colCount, stdDevs, colCount,
			modelParams, modelParamsCount,
			precomputedWeights, precomputedWeightsCount)
		<< endl;

	cout << "\nPredicting: "
		<< predictByLSSVM(kernel, kernelParams, kernelParamsCount,
			normalizedDataset, colCount, rowCount,
			labels, labelsLength,
			means, colCount,
			stdDevs, colCount,
			modelParams, modelParamsCount,
			precomputedWeights, precomputedWeightsCount,
			dataset, rowCount, colCount,
			predictedLabels, labelsLength)
		<< endl;

	// Saving

	cout << "\nSaving data ...\n";

	string fileName = "dataset.csv";

	ofstream file(fileName, ios::out);
	if (!file)
	{
		cout << "FAIL TO CREATE FILE FOR RESULTS!\n";
		return;
	}

	for (int j = 1; j <= colCount; j++)
		file << 'x' << j << ',';
	file << "y,predicted,faults\n";

	for (int i = 0; i < rowCount; i++)
	{
		for (int j = 0; j < colCount; j++)
			file << dataset[i + j * rowCount] << ',';

		file << labels[i] << ","
			<< predictedLabels[i] << ","
			<< labels[i] * predictedLabels[i] << endl;
	}

	delete[] kernelParams;
	delete[] dataset;
	delete[] labels;
	delete[] normalizedDataset;
	delete[] means;
	delete[] stdDevs;
	delete[] modelParams;
	delete[] precomputedWeights;
	delete[] predictedLabels;
}
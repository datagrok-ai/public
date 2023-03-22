// tests.cpp

// Implementations of funtions for testing methods.

#include<iostream>
#include<cstring>
#include<fstream>
#include <chrono>
using namespace std;
using namespace std::chrono;

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

void testCreateNormalizedDataset()
{
	const int ROW_COUNT = 5;
	const int COL_COUNT = 3;
	const int SIZE = ROW_COUNT * COL_COUNT;

	float data[SIZE] = { 0 };
	float means[COL_COUNT] = { 0 };
	float stdDevs[COL_COUNT] = { 0 };

	int index = 0;

	for (int j = 0; j < COL_COUNT; j++)
	{
		for (int i = 0; i < ROW_COUNT; i++)
			data[index++] = (float)(2 * j + (j + 1) * i);
	}

	cout << "Data:\n";
	for (int i = 0; i < SIZE; i++)
		cout << "  " << data[i];

	float dataset[SIZE] = { 0 };

	cout << "\nDataset before:\n";
	for (int i = 0; i < SIZE; i++)
		cout << "  " << dataset[i];

	createNormalizedDataset(data, ROW_COUNT, COL_COUNT, dataset, means, stdDevs);

	cout << "\nDataset after:\n";
	for (int i = 0; i < SIZE; i++)
		cout << "  " << dataset[i];

	cout << "\nMeans:\n";
	for (int i = 0; i < COL_COUNT; i++)
		cout << "  " << means[i];

	cout << "\nStandard deviations:\n";
	for (int i = 0; i < COL_COUNT; i++)
		cout << "  " << stdDevs[i];

} // testCreateNormalizedDataset

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

void testTrainModelNormalizedDataLinear()
{
	const int featuresCount = 2;
	const int samplesCount = 10000;

	float w[featuresCount] = { 1, -1 };
	float b = 2;

	/*float minVal[featuresCount] = { -150, -2700 };
	float maxVal[featuresCount] = { 2500, 1700 };*/

	float minVal[featuresCount] = { -500, -600};
	float maxVal[featuresCount] = { 500, 600};

	/*float minVal[featuresCount] = {-1, -1};
	float maxVal[featuresCount] = { 1, 1 };*/

	float violatorsPercentage = 30;
	float actualViolatorsPercentage = 0;

	string fileName = "svm_test.csv";

	float* data = new float[featuresCount * samplesCount];
	float means[featuresCount] = { 0 };
	float stdDevs[featuresCount] = { 0 };
	float* labels = new float[samplesCount];
	float* xTrain = new float[featuresCount * samplesCount];
	float* yTrain = labels;

	/*cout << "Generate data: "
		<< generateLinearSeparable(w, b, featuresCount, samplesCount,
			minVal, maxVal, data, labels) << endl;*/

	cout << "Generate data: "
		<< generateLinearNonSeparable(w, b, featuresCount, samplesCount,
			minVal, maxVal, data, labels, violatorsPercentage, actualViolatorsPercentage) << endl;

	cout << "\nActual violators, %: " << actualViolatorsPercentage << endl;

	//cout << "\nCreate xTrain: " << createDataset(data, samplesCount, featuresCount, xTrain) << endl;

	cout << "\nCreate xTrain: " << createNormalizedDataset(data, samplesCount, featuresCount, xTrain, means, stdDevs) << endl;

	/*cout << "\nX train:\n";
	for (int i = 0; i < featuresCount * samplesCount; i++)
		cout << "  " << xTrain[i];

	cout << "\nyTrain:\n";
	for (int i = 0; i < samplesCount; i++)
		cout << "  " << yTrain[i];*/

	float gamma = 2.0f;
	int kernel = LINEAR;
	float kernelParams[MAX_NUM_OF_KERNEL_PARAM] = { 1.0f, 1.0f };

	float modelParams[samplesCount + 1] = { 0.0f };

	cout << "\n\nTrain model result: "
		<< trainLSSVM(gamma, kernel, kernelParams, xTrain, yTrain, samplesCount, featuresCount, modelParams)
		<< endl;

	/*cout << "\nParameters are:\n";
	for (int i = 0; i <= samplesCount; i++)
		cout << "  " << modelParams[i];

	cout << "\nbias: " << modelParams[samplesCount] << endl;*/

	// PREDICTION
	
	// simple
	/*const int targetSamplesCount = 4;
	float targetData[featuresCount * targetSamplesCount] = {-200, 200, 100, -50,
		200, -200, 120, -170};
	float targetLabels[targetSamplesCount] = { 0 };

	cout << "\nPrediction: "
		<< predictByLSSVM(kernel, kernelParams,
			xTrain, yTrain, samplesCount, featuresCount,
			means, stdDevs, modelParams,
			targetData, targetLabels, targetSamplesCount) << endl;*/

	float* predictedLabels = new float[samplesCount];
	float error = 0;

	cout << "\nPredicting: "
		<< predictByLSSVM(kernel, kernelParams,
			xTrain, yTrain, samplesCount, featuresCount,
			means, stdDevs, modelParams,
			data, predictedLabels, samplesCount) << endl;

	/*cout << "\nPrediction:\n";

	for (int i = 0; i < samplesCount; i++)
		cout << "  " << yTrain[i] << "  <->  " << predictedLabels[i] << endl;*/

	evalutaeErrorPercent(yTrain, predictedLabels, samplesCount, error);

	cout << "-----------------------------------------\n"
		<< "Error, %: " << error << endl;


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
	delete[] predictedLabels;
} // testTrainModelNormalizedDataLinear

void testTrainModelNormalizedDataLinearHighDim(bool saveData)
{
	const int featuresCount = 2;
	const int samplesCount = 20000;
	const int samplesCountTest = samplesCount / 10;

	float* w = new float [featuresCount];
	float b = 0;

	cout << "Generating model parameters: "
		<< generateModelParams(w, b, featuresCount) << endl;

	/*cout << "\nw:\n";
	for (int i = 0; i < featuresCount; i++)
		cout << " " << w[i];

	cout << "\nbias: " << b << endl;*/

	/*float minVal[featuresCount] = { -150, -2700 };
	float maxVal[featuresCount] = { 2500, 1700 };*/

	float min = -450;
	float max = 560;

	float* minVal = new float [featuresCount];
	float* maxVal = new float [featuresCount];

	for (int i = 0; i < featuresCount; i++)
	{
		minVal[i] = min;
		maxVal[i] = max;
	}

	/*float minVal[featuresCount] = {-1, -1};
	float maxVal[featuresCount] = { 1, 1 };*/

	float violatorsPercentage = 0;
	float actualViolatorsPercentageTrain = 0;
	float actualViolatorsPercentageTest = 0;

	string fileNameTrain = "svm_train.csv";
	string fileNameTest = "svm_test.csv";

	float* dataTrain = new float[featuresCount * samplesCount];
	float* dataTest = new float[featuresCount * samplesCountTest];
	float* means = new float [featuresCount];
	float* stdDevs = new float [featuresCount];
	float* labels = new float[samplesCount];
	float* labelsTest = new float[samplesCountTest];
	float* xTrain = new float[featuresCount * samplesCount];
	float* yTrain = labels;

	/*cout << "Generate data: "
		<< generateLinearSeparable(w, b, featuresCount, samplesCount,
			minVal, maxVal, data, labels) << endl;*/

	auto start = high_resolution_clock::now();

	cout << "\nGenerate train data: "
		<< generateLinearNonSeparable(w, b, featuresCount, samplesCount,
			minVal, maxVal, dataTrain, labels, violatorsPercentage, actualViolatorsPercentageTrain) << endl;

	auto finish = high_resolution_clock::now();

	auto duration = duration_cast<microseconds>(finish - start);

	cout << "Time: " << 1e-6 * duration.count() << " sec.\n";

	cout << "\nActual train violators, %: " << actualViolatorsPercentageTrain << endl;

	cout << "\nGenerate test data: "
		<< generateLinearNonSeparable(w, b, featuresCount, samplesCountTest,
			minVal, maxVal, dataTest, labelsTest, violatorsPercentage, actualViolatorsPercentageTest) << endl;

	cout << "\nActual test violators, %: " << actualViolatorsPercentageTest << endl;

	//cout << "\nCreate xTrain: " << createDataset(data, samplesCount, featuresCount, xTrain) << endl;

	start = high_resolution_clock::now();

	cout << "\nCreate xTrain: " << createNormalizedDataset(dataTrain, samplesCount, featuresCount, xTrain, means, stdDevs) << endl;

	finish = high_resolution_clock::now();

	auto durationData = duration_cast<microseconds>(finish - start);

	cout << "Time: " << 1e-6 * durationData.count() << " sec.\n";

	/*cout << "\nX train:\n";
	for (int i = 0; i < featuresCount * samplesCount; i++)
		cout << "  " << xTrain[i];

	cout << "\nyTrain:\n";
	for (int i = 0; i < samplesCount; i++)
		cout << "  " << yTrain[i];*/

	float gamma = 1.0f;
	int kernel = LINEAR;
	float kernelParams[MAX_NUM_OF_KERNEL_PARAM] = { 1.0f, 1.0f };

	float* modelParams = new float [samplesCount + 1];

	start = high_resolution_clock::now();

	cout << "\nTrain model result: "
		<< trainLSSVM(gamma, kernel, kernelParams, xTrain, yTrain, samplesCount, featuresCount, modelParams)
		<< endl;

	finish = high_resolution_clock::now();

	auto durationTrain = duration_cast<microseconds>(finish - start);

	cout << "Time: " << 1e-6 * durationTrain.count() << " sec.\n";
	cout << "TRAIN MODEL TIME: " << 1e-6 * durationTrain.count() + 1e-6 * durationData.count() << " sec.\n";

	/*cout << "\nParameters are:\n";
	for (int i = 0; i <= samplesCount; i++)
		cout << "  " << modelParams[i];

	cout << "\nbias: " << modelParams[samplesCount] << endl;*/

	// PREDICTION	

	float* predictedLabelsTrain = new float[samplesCount];
	float errorTrain = 0;

	start = high_resolution_clock::now();

	cout << "\nPredicting train: "
		<< predictByLSSVM(kernel, kernelParams,
			xTrain, yTrain, samplesCount, featuresCount,
			means, stdDevs, modelParams,
			dataTrain, predictedLabelsTrain, samplesCount) << endl;

	finish = high_resolution_clock::now();

	duration = duration_cast<microseconds>(finish - start);

	cout << "Time: " << 1e-6 * duration.count() << " sec.\n";

	/*cout << "\nPrediction:\n";

	for (int i = 0; i < samplesCount; i++)
		cout << "  " << yTrain[i] << "  <->  " << predictedLabels[i] << endl;*/

	evalutaeErrorPercent(yTrain, predictedLabelsTrain, samplesCount, errorTrain);

	float* predictedLabelsTest = new float[samplesCountTest];
	float errorTest = 0;

	cout << "\nPredicting test: "
		<< predictByLSSVM(kernel, kernelParams,
			xTrain, yTrain, samplesCount, featuresCount,
			means, stdDevs, modelParams,
			dataTest, predictedLabelsTest, samplesCountTest) << endl;

	/*cout << "\nPrediction:\n";

	for (int i = 0; i < samplesCount; i++)
		cout << "  " << yTrain[i] << "  <->  " << predictedLabels[i] << endl;*/

	evalutaeErrorPercent(labelsTest, predictedLabelsTest, samplesCountTest, errorTest);

	cout << "-----------------------------------------\n"
		<< "Train error, %: " << errorTrain << endl
		<< "Test error, %: " << errorTest << endl;
	
	if (!saveData)
		return;

	cout << "\nSaving to file ...\n";

	ofstream fileTrain(fileNameTrain, ios::out);
	if (!fileTrain)
	{
		cout << "FAIL TO CREATE FILE FOR RESULTS!\n";
		return;
	}

	for (int j = 1; j <= featuresCount; j++)
		fileTrain << 'x' << j << ',';
	fileTrain << "label,predicted,fault\n";

	for (int i = 0; i < samplesCount; i++)
	{
		for (int j = 0; j < featuresCount; j++)
			fileTrain << dataTrain[i + j * samplesCount] << ',';

		fileTrain << labels[i] << ","
			<< predictedLabelsTrain[i] << ","
			<< labels[i] * predictedLabelsTrain[i] << endl;
	}

	fileTrain.close();

	ofstream fileTest(fileNameTest, ios::out);
	if (!fileTest)
	{
		cout << "FAIL TO CREATE FILE FOR RESULTS!\n";
		return;
	}

	for (int j = 1; j <= featuresCount; j++)
		fileTest << 'x' << j << ',';
	fileTest << "label,predicted,fault\n";

	for (int i = 0; i < samplesCountTest; i++)
	{
		for (int j = 0; j < featuresCount; j++)
			fileTest << dataTest[i + j * samplesCountTest] << ',';

		fileTest << labelsTest[i] << ","
			<< predictedLabelsTest[i] << ","
			<< labelsTest[i] * predictedLabelsTest[i] << endl;
	}

	fileTest.close();

	delete[] dataTrain;
	delete[] dataTest;
	delete[] labels;
	delete[] labelsTest;
	delete[] xTrain;
	delete[] predictedLabelsTrain;
	delete[] predictedLabelsTest;
	delete[] minVal;
	delete[] maxVal;
	delete[] w;
	delete[] means;
	delete[] stdDevs;
	delete[] modelParams;
} // testTrainModelNormalizedDataLinearHighDim

void testTrainModelNormalizedDataLinearHighDimDouble()
{
	const int featuresCount = 10000;
	const int samplesCount = 10000;
	const int samplesCountTest = samplesCount / 10;

	double* w = new double[featuresCount];
	double b = 0;

	cout << "Generating model parameters: "
		<< generateModelParams(w, b, featuresCount) << endl;

	/*cout << "\nw:\n";
	for (int i = 0; i < featuresCount; i++)
		cout << " " << w[i];

	cout << "\nbias: " << b << endl;*/

	/*float minVal[featuresCount] = { -150, -2700 };
	float maxVal[featuresCount] = { 2500, 1700 };*/

	double min = -450;
	double max = 560;

	double* minVal = new double[featuresCount];
	double* maxVal = new double[featuresCount];

	for (int i = 0; i < featuresCount; i++)
	{
		minVal[i] = min;
		maxVal[i] = max;
	}

	/*float minVal[featuresCount] = {-1, -1};
	float maxVal[featuresCount] = { 1, 1 };*/

	double violatorsPercentage = 0;
	double actualViolatorsPercentageTrain = 0;
	double actualViolatorsPercentageTest = 0;

	string fileName = "svm_train.csv";

	double* dataTrain = new double[featuresCount * samplesCount];
	double* dataTest = new double[featuresCount * samplesCountTest];
	double* means = new double[featuresCount];
	double* stdDevs = new double[featuresCount];
	double* labels = new double[samplesCount];
	double* labelsTest = new double[samplesCountTest];
	double* xTrain = new double[featuresCount * samplesCount];
	double* yTrain = labels;

	/*cout << "Generate data: "
		<< generateLinearSeparable(w, b, featuresCount, samplesCount,
			minVal, maxVal, data, labels) << endl;*/

	cout << "\nGenerate train data: "
		<< generateLinearNonSeparable(w, b, featuresCount, samplesCount,
			minVal, maxVal, dataTrain, labels, violatorsPercentage, actualViolatorsPercentageTrain) << endl;

	cout << "\nActual train violators, %: " << actualViolatorsPercentageTrain << endl;

	cout << "\nGenerate test data: "
		<< generateLinearNonSeparable(w, b, featuresCount, samplesCountTest,
			minVal, maxVal, dataTest, labelsTest, violatorsPercentage, actualViolatorsPercentageTest) << endl;

	cout << "\nActual test violators, %: " << actualViolatorsPercentageTest << endl;

	//cout << "\nCreate xTrain: " << createDataset(data, samplesCount, featuresCount, xTrain) << endl;

	cout << "\nCreate xTrain: " << createNormalizedDataset(dataTrain, samplesCount, featuresCount, xTrain, means, stdDevs) << endl;

	/*cout << "\nX train:\n";
	for (int i = 0; i < featuresCount * samplesCount; i++)
		cout << "  " << xTrain[i];

	cout << "\nyTrain:\n";
	for (int i = 0; i < samplesCount; i++)
		cout << "  " << yTrain[i];*/

	double gamma = 1.0f;
	int kernel = LINEAR;
	double kernelParams[MAX_NUM_OF_KERNEL_PARAM] = { 1.0, 1.0 };

	double* modelParams = new double[samplesCount + 1];

	cout << "\n\nTrain model result: "
		<< trainLSSVM(gamma, kernel, kernelParams, xTrain, yTrain, samplesCount, featuresCount, modelParams)
		<< endl;

	/*cout << "\nParameters are:\n";
	for (int i = 0; i <= samplesCount; i++)
		cout << "  " << modelParams[i];

	cout << "\nbias: " << modelParams[samplesCount] << endl;*/

	// PREDICTION	

	double* predictedLabelsTrain = new double[samplesCount];
	double errorTrain = 0;

	cout << "\nPredicting train: "
		<< predictByLSSVM(kernel, kernelParams,
			xTrain, yTrain, samplesCount, featuresCount,
			means, stdDevs, modelParams,
			dataTrain, predictedLabelsTrain, samplesCount) << endl;

	/*cout << "\nPrediction:\n";

	for (int i = 0; i < samplesCount; i++)
		cout << "  " << yTrain[i] << "  <->  " << predictedLabels[i] << endl;*/

	evalutaeErrorPercent(yTrain, predictedLabelsTrain, samplesCount, errorTrain);


	double* predictedLabelsTest = new double[samplesCountTest];
	double errorTest = 0;

	cout << "\nPredicting test: "
		<< predictByLSSVM(kernel, kernelParams,
			xTrain, yTrain, samplesCount, featuresCount,
			means, stdDevs, modelParams,
			dataTest, predictedLabelsTest, samplesCountTest) << endl;

	/*cout << "\nPrediction:\n";

	for (int i = 0; i < samplesCount; i++)
		cout << "  " << yTrain[i] << "  <->  " << predictedLabels[i] << endl;*/

	evalutaeErrorPercent(labelsTest, predictedLabelsTest, samplesCountTest, errorTest);

	cout << "-----------------------------------------\n"
		<< "Train error, %: " << errorTrain << endl
		<< "Test error, %: " << errorTest << endl;

	cout << "\nSaving to file ...\n";

	ofstream file(fileName, ios::out);
	if (!file)
	{
		cout << "FAIL TO CREATE FILE FOR RESULTS!\n";
		return;
	}

	for (int j = 1; j <= featuresCount; j++)
		file << 'x' << j << ',';
	file << "label,predicted,fault\n";

	for (int i = 0; i < samplesCount; i++)
	{
		for (int j = 0; j < featuresCount; j++)
			file << dataTrain[i + j * samplesCount] << ',';

		file << labels[i] << ","
			<< predictedLabelsTrain[i] << ","
			<< labels[i] * predictedLabelsTrain[i] << endl;
	}

	delete[] dataTrain;
	delete[] dataTest;
	delete[] labels;
	delete[] labelsTest;
	delete[] xTrain;
	delete[] predictedLabelsTrain;
	delete[] predictedLabelsTest;
	delete[] minVal;
	delete[] maxVal;
	delete[] w;
	delete[] means;
	delete[] stdDevs;
	delete[] modelParams;
} // testTrainModelNormalizedDataLinearHighDimDouble
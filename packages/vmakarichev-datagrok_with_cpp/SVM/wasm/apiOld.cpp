// api.cpp

// Implementation of API-functions for DATAGROK packages.

#include <emscripten.h>

// The following provides convenient naming of the exported functions.
extern "C" {
	
    int generateDataset(int featuresCount, int samplesCount,	
	    float min, float max, float violatorsPercentage,
	    float * dataset, int rowCount, int colCount,
	    float * labels, int labelsLength);

	int demo(float gamma,
		float * datasetTrain, int datasetTrainRowCount, int datasetTrainColCount,
	    float * labelsTrain, int labelsTrainLength,
	    float * labelsTrainPredicted, int labelsTrainPredictedLength);
	
	int demoUPD(float gamma,
		float * datasetTrain, int datasetTrainRowCount, int datasetTrainColCount,
		float * labelsTrain, int labelsTrainLength,
		float * datasetTest, int datasetTestRowCount, int datasetTestColCount,
		float * labelsTrainPredicted, int labelsTrainPredictedLength,
		float * labelsTestPredicted, int labelsTestPredictedLength);

	int demoUPDupdUPD(float gamma,
		float * datasetTrain, int datasetTrainRowCount, int datasetTrainColCount,
		float * labelsTrain, int labelsTrainLength,
		float * datasetTest, int datasetTestRowCount, int datasetTestColCount,
		float * labelsTrainPredicted, int labelsTrainPredictedLength,
		float * labelsTestPredicted, int labelsTestPredictedLength);

	int demoLinear(float gamma,
		float * datasetTrain, int datasetTrainRowCount, int datasetTrainColCount,
		float * labelsTrain, int labelsTrainLength,
		float * datasetTest, int datasetTestRowCount, int datasetTestColCount,
		int weightsCount,
		int paramsCount,
		float * labelsTrainPredicted, int labelsTrainPredictedLength,
		float * labelsTestPredicted, int labelsTestPredictedLength,
		float * means, int meansLength,
		float * stdDevs, int stdDevsLength,
		float * weights, int weightsLength,
		float * modelParams, int modelParamsLength);
}

#include "generators.h"
#include "svm.h"

//name: generateDataset
//input: int featuresCount
//input: int samplesCount
//input: double min 
//input: double max
//input: double violatorsPercentage
//output: column_list dataset [new(samplesCount, featuresCount)]
//output: column labels [new(samplesCount)]
EMSCRIPTEN_KEEPALIVE
int generateDataset(int featuresCount, int samplesCount,	
	float min, float max, float violatorsPercentage,
	float * dataset, int rowCount, int colCount,
	float * labels, int labelsLength)
{
	using namespace gener;

	float* w = new float[featuresCount];
	float b = 0;

	int resCode = generateModelParams(w, b, featuresCount);
	if (resCode != NO_ERRORS)
		return resCode;

	float* minVal = new float[featuresCount];
	float* maxVal = new float[featuresCount];

	for (int i = 0; i < featuresCount; i++)
	{
		minVal[i] = min;
		maxVal[i] = max;
	}

	float actualViolatorsPercentage = 0;

	return generateLinearNonSeparable(w, b, featuresCount, samplesCount,
			minVal, maxVal, dataset, labels, violatorsPercentage, actualViolatorsPercentage);
} // generateDataset

//name: demoSVM
//input: double gamma
//input: dataframe dfTrain
//input: column_list datasetTrain
//input: column labelsTrain
//output: column labelsTrainPredicted [new(labelsTrain.rowCount)]
EMSCRIPTEN_KEEPALIVE
int demo(float gamma,
	float * datasetTrain, int datasetTrainRowCount, int datasetTrainColCount,
	float * labelsTrain, int labelsTrainLength,
	float * labelsTrainPredicted, int labelsTrainPredictedLength)
{
	using namespace svm;

	//float gamma = 1.0f;

	int kernel = LINEAR;
	float kernelParams[MAX_NUM_OF_KERNEL_PARAM] = { 1.0f, 1.0f };
	int samplesCount = datasetTrainRowCount;
	int featuresCount = datasetTrainColCount;

	float* xTrain = new float[featuresCount * samplesCount];
	float* yTrain = labelsTrain;
	float* means = new float[featuresCount];
	float* stdDevs = new float[featuresCount];
	float* modelParams = new float[samplesCount + 1];
		
	createNormalizedDataset(datasetTrain, samplesCount, featuresCount, xTrain, means, stdDevs);
	
	trainLSSVM(gamma, kernel, kernelParams, xTrain, yTrain, samplesCount, featuresCount, modelParams);	
	
	predictByLSSVM(kernel, kernelParams,
		xTrain, yTrain, samplesCount, featuresCount,
		means, stdDevs, modelParams,
		datasetTrain, labelsTrainPredicted, samplesCount);

	delete[] xTrain;
	delete[] means;
	delete[] stdDevs;
	delete[] modelParams;
	
	return 0;
} // demo

//name: demoSVMupd
//input: double gamma
//input: dataframe dfTrain
//input: column_list datasetTrain
//input: column labelsTrain
//input: dataframe dfTest
//input: column_list datasetTest
//output: column labelsTrainPredicted [new(datasetTrain.rowCount)]
//output: column labelsTestPredicted [new(datasetTest.rowCount)]
EMSCRIPTEN_KEEPALIVE
int demoUPD(float gamma,
	float * datasetTrain, int datasetTrainRowCount, int datasetTrainColCount,
	float * labelsTrain, int labelsTrainLength,
	float * datasetTest, int datasetTestRowCount, int datasetTestColCount,
	float * labelsTrainPredicted, int labelsTrainPredictedLength,
	float * labelsTestPredicted, int labelsTestPredictedLength)
{
	using namespace svm;

	int kernel = LINEAR;
	float kernelParams[MAX_NUM_OF_KERNEL_PARAM] = { 1.0f, 1.0f };
	int samplesCount = datasetTrainRowCount;
	int samplesCountTest = datasetTestRowCount;
	int featuresCount = datasetTrainColCount;

	float* xTrain = new float[featuresCount * samplesCount];
	float* yTrain = labelsTrain;
	float* means = new float[featuresCount];
	float* stdDevs = new float[featuresCount];
	float* modelParams = new float[samplesCount + 1];
		
	createNormalizedDataset(datasetTrain, samplesCount, featuresCount, xTrain, means, stdDevs);
	
	trainLSSVM(gamma, kernel, kernelParams, xTrain, yTrain, samplesCount, featuresCount, modelParams);	
	
	predictByLSSVM(kernel, kernelParams,
		xTrain, yTrain, samplesCount, featuresCount,
		means, stdDevs, modelParams,
		datasetTrain, labelsTrainPredicted, samplesCount);

	predictByLSSVM(kernel, kernelParams,
		xTrain, yTrain, samplesCount, featuresCount,
		means, stdDevs, modelParams,
		datasetTest, labelsTestPredicted, samplesCountTest);

	delete[] xTrain;
	delete[] means;
	delete[] stdDevs;
	delete[] modelParams;
	
	return 0;
} // demoUPD

//name: demoSVMupdUPD
//input: double gamma
//input: dataframe dfTrain
//input: column_list datasetTrain
//input: column labelsTrain
//input: dataframe dfTest
//input: column_list datasetTest
//output: column labelsTrainPredicted [new(datasetTrain.rowCount)]
//output: column labelsTestPredicted [new(datasetTest.rowCount)]
EMSCRIPTEN_KEEPALIVE
int demoUPDupdUPD(float gamma,
	float * datasetTrain, int datasetTrainRowCount, int datasetTrainColCount,
	float * labelsTrain, int labelsTrainLength,
	float * datasetTest, int datasetTestRowCount, int datasetTestColCount,
	float * labelsTrainPredicted, int labelsTrainPredictedLength,
	float * labelsTestPredicted, int labelsTestPredictedLength)
{
	using namespace svm;

	int kernel = LINEAR;
	float kernelParams[MAX_NUM_OF_KERNEL_PARAM] = { 1.0f, 1.0f };
	int samplesCount = datasetTrainRowCount;
	int samplesCountTest = datasetTestRowCount;
	int featuresCount = datasetTrainColCount;

	float* xTrain = new float[featuresCount * samplesCount];
	float* yTrain = labelsTrain;
	float* means = new float[featuresCount];
	float* stdDevs = new float[featuresCount];
	float* modelParams = new float[samplesCount + 1];
	float* weights = new float[featuresCount + 1];
		
	createNormalizedDataset(datasetTrain, samplesCount, featuresCount, xTrain, means, stdDevs);
	
	trainLSSVM(gamma, kernel, kernelParams, 
	    xTrain, yTrain, samplesCount, featuresCount, 
		modelParams, weights);	
	
	predictByLSSVM(kernel, kernelParams,
		xTrain, yTrain, samplesCount, featuresCount,
		means, stdDevs, modelParams, weights,
		datasetTrain, labelsTrainPredicted, samplesCount);

	predictByLSSVM(kernel, kernelParams,
		xTrain, yTrain, samplesCount, featuresCount,
		means, stdDevs, modelParams, weights,
		datasetTest, labelsTestPredicted, samplesCountTest);

	delete[] xTrain;
	delete[] means;
	delete[] stdDevs;
	delete[] modelParams;
	delete[] weights;
	
	return 0;
} // demoUPDupd

//name: demoLinear
//input: double gamma
//input: dataframe dfTrain
//input: column_list datasetTrain
//input: column labelsTrain
//input: dataframe dfTest
//input: column_list datasetTest
//input: int weightsCount
//input: int paramsCount
//output: column labelsTrainPredicted [new(datasetTrain.rowCount)]
//output: column labelsTestPredicted [new(datasetTest.rowCount)]
//output: column means [new(datasetTrain.columnCount)]
//output: column stdDevs [new(datasetTrain.columnCount)]
//output: column weights [new(weightsCount)]
//output: column modelParams [new(paramsCount)]
EMSCRIPTEN_KEEPALIVE
int demoLinear(float gamma,
	float * datasetTrain, int datasetTrainRowCount, int datasetTrainColCount,
	float * labelsTrain, int labelsTrainLength,
	float * datasetTest, int datasetTestRowCount, int datasetTestColCount,
	int weightsCount,
	int paramsCount,
	float * labelsTrainPredicted, int labelsTrainPredictedLength,
	float * labelsTestPredicted, int labelsTestPredictedLength,
	float * means, int meansLength,
	float * stdDevs, int stdDevsLength,
	float * weights, int weightsLength,
	float * modelParams, int modelParamsLength)
{
	using namespace svm;

	int kernel = LINEAR;
	float kernelParams[MAX_NUM_OF_KERNEL_PARAM] = { 1.0f, 1.0f };
	int samplesCount = datasetTrainRowCount;
	int samplesCountTest = datasetTestRowCount;
	int featuresCount = datasetTrainColCount;

	float* xTrain = new float[featuresCount * samplesCount];
	float* yTrain = labelsTrain;
		
	createNormalizedDataset(datasetTrain, samplesCount, featuresCount, xTrain, means, stdDevs);
	
	trainLSSVM(gamma, kernel, kernelParams, 
	    xTrain, yTrain, samplesCount, featuresCount, 
		modelParams, weights);	
	
	predictByLSSVM(kernel, kernelParams,
		xTrain, yTrain, samplesCount, featuresCount,
		means, stdDevs, modelParams, weights,
		datasetTrain, labelsTrainPredicted, samplesCount);

	predictByLSSVM(kernel, kernelParams,
		xTrain, yTrain, samplesCount, featuresCount,
		means, stdDevs, modelParams, weights,
		datasetTest, labelsTestPredicted, samplesCountTest);

	delete[] xTrain;
	
	return 0;
} // demoLinear

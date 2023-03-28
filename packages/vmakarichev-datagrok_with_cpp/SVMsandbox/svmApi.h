// svmApi.h

#ifndef SVM_API_H
#define SVM_API_H

#include "svm.h"

// API-functions

/*
*/
int generateDataset(int kernel, 
	float * kernelParams, int kernelParamsCount,
	float * dataset, int rowCount, int colCount,
	float * labels, int labelsLength,
	float min, float max, float violatorsPercentage);

int normalizeDataset(float * dataset, int datasetRowCount, int datasetColCount,
	float * normalizedData, int normalizedDataRowCount, int normalizedDataColCount,
	float * means, int meansLength,
	float * stdDevs, int stdDevsLength);

int trainLSSVM(float gamma, int kernel,
	float * kernelParams, int kernelParamsCount,
	int modelParamsCount, int precomputedWeightsCount,
	float * dataset, int datasetRowCount, int datasetColCount,
	float * labels, int labelsLength,
	float * normalizedData, int normalizedDataRowCount, int normalizedDataColCount,
	float * means, int meansLength,
	float * stdDevs, int stdDevsLength,
	float * modelParams, int modelParamsLength,
	float * precomputedWeights, int precomputedWeightsLength);

int predictByLSSVM(int kernel,
	float* kernelParams, int kernelParamsCount,	
	float* normalizedData, int normalizedDataRowCount, int normalizedDataColCount,
	float* labels, int labelsLength,
	float* means, int meansLength,
	float* stdDevs, int stdDevsLength,
	float* modelParams, int modelParamsLength,
	float* precomputedWeights, int precomputedWeightsLength,
	float* targetData, int targetDataRowCount, int targetDataColCount,
	float* prediction, int predictionLength);


#endif // SVM_API_H


#include "svm.h"
#include "dataGenerators.h"
#include "dataMining.h"
#include "svmApi.h"

int generateDataset(int kernel,
	float* kernelParams, int kernelParamsCount,
	float* dataset, int rowCount, int colCount,
	float* labels, int labelsLength,
	float min, float max, float violatorsPercentage)
{
	using namespace svm;

	if (kernelParamsCount != MAX_NUM_OF_KERNEL_PARAM)
		return INCORRECT_KERNEL_PARAMS_COUNT;

	if (!areKernelParametersCorrect(kernel, kernelParams))
		return INCORRECT_PARAMETER_OF_KERNEL;

	if ((rowCount < 1) || (colCount < 1) || (rowCount != labelsLength))
		return INCORRECT_SIZE;

	return generateNonSeparable(kernel, kernelParams,
		colCount, rowCount, min, max, 
		dataset, labels, violatorsPercentage);
} // generateDataset

int normalizeDataset(float * dataset, int datasetRowCount, int datasetColCount,
	float * normalizedData, int normalizedDataRowCount, int normalizedDataColCount,
	float * means, int meansLength,
	float * stdDevs, int stdDevsLength)
{
	using namespace dmt;

	// check sizes
	if ((datasetRowCount != normalizedDataColCount)
		|| (datasetColCount != normalizedDataRowCount)
		|| (meansLength != datasetColCount)
		|| (stdDevsLength != datasetColCount))
		return INCORRECT_SIZE;

	return getNormalizedDataset(dataset, datasetRowCount, datasetColCount,
		normalizedData, means, stdDevs);
}

int trainLSSVM(float gamma, int kernel,
	float * kernelParams, int kernelParamsCount,
	int modelParamsCount, int precomputedWeightsCount,
	float * dataset, int datasetRowCount, int datasetColCount,
	float * labels, int labelsLength,
	float * normalizedData, int normalizedDataRowCount, int normalizedDataColCount,
	float * means, int meansLength,
	float * stdDevs, int stdDevsLength,
	float * modelParams, int modelParamsLength,
	float * precomputedWeights, int precomputedWeightsLength)
{
	using namespace svm;
	using dmt::getNormalizedDataset;

	// check gamma 
	if (!isGammaCorrect(gamma))
		return INCORRECT_HYPERPARAMETER;

	// check kernel params count
	if (kernelParamsCount != MAX_NUM_OF_KERNEL_PARAM)
		return INCORRECT_KERNEL_PARAMS_COUNT;

	// check kernel specification
	if (!areKernelParametersCorrect(kernel, kernelParams))
		return INCORRECT_PARAMETER_OF_KERNEL;

	// check sizes
	if ((datasetRowCount < 1)
		|| (datasetColCount < 1)
		|| (labelsLength != datasetRowCount)
		|| (normalizedDataRowCount != datasetColCount)
		|| (normalizedDataColCount != datasetRowCount)
		|| (meansLength != datasetColCount)
		|| (stdDevsLength != datasetColCount)
		|| (modelParamsLength != datasetRowCount + 1)
		|| (precomputedWeightsLength != datasetColCount + 1)) 
		return INCORRECT_SIZE;

	// normalize data
	int resCode = getNormalizedDataset(dataset, datasetRowCount, datasetColCount,
		normalizedData, means, stdDevs);
	if (resCode != dmt::NO_ERRORS)
		return resCode;

	// train LS-SVM model
	return trainLSSVM(gamma, kernel, kernelParams,
		normalizedData, labels, datasetRowCount, datasetColCount,
		modelParams, precomputedWeights);
} // trainLSSVM

int predictByLSSVM(int kernel,
	float* kernelParams, int kernelParamsCount,
	float* normalizedData, int normalizedDataRowCount, int normalizedDataColCount,
	float* labels, int labelsLength,
	float* means, int meansLength,
	float* stdDevs, int stdDevsLength,
	float* modelParams, int modelParamsLength,
	float* precomputedWeights, int precomputedWeightsLength,
	float* targetData, int targetDataRowCount, int targetDataColCount,
	float* prediction, int predictionLength)
{
	using namespace svm;

	// check kernel params count
	if (kernelParamsCount != MAX_NUM_OF_KERNEL_PARAM)
		return INCORRECT_KERNEL_PARAMS_COUNT;

	// check kernel specification
	if (!areKernelParametersCorrect(kernel, kernelParams))
		return INCORRECT_PARAMETER_OF_KERNEL;
	
	// check sizes
	if ((normalizedDataRowCount < 1)
		|| (normalizedDataColCount < 1)
		|| (labelsLength != normalizedDataColCount) 
		|| (meansLength != normalizedDataRowCount)
		|| (stdDevsLength != normalizedDataRowCount)
		|| (modelParamsLength != normalizedDataColCount + 1)
		|| (precomputedWeightsLength != normalizedDataRowCount + 1)
		|| (targetDataRowCount < 1)
		|| (targetDataColCount < 1)
		|| (targetDataColCount != normalizedDataRowCount)
		|| (predictionLength != targetDataRowCount))
		return INCORRECT_SIZE;

	// predict labels
	return predictByLSSVM(kernel, kernelParams, normalizedData, labels,
		normalizedDataColCount, normalizedDataRowCount,
		means, stdDevs, modelParams, precomputedWeights,
		targetData, prediction, targetDataRowCount);
} // predictByLSSVM
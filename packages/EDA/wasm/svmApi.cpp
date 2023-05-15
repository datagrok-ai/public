#include "svm.h"
#include "dataGenerators.h"
#include "dataMining.h"

#include <emscripten.h>

// The following provides convenient naming of the exported functions.
extern "C" {
	int generateDataset(int kernel,
		float * kernelParams, int kernelParamsCount,
		int samplesCount, int featuresCount,
		float min, float max,
		float violatorsPercentage,
		float * dataset, int rowCount, int colCount,
		float * labels, int labelsLength) noexcept;		

	int normalizeDataset(float * dataset, int datasetRowCount, int datasetColCount,
		float * normalizedData, int normalizedDataRowCount, int normalizedDataColCount,
		float * means, int meansLength,
		float * stdDevs, int stdDevsLength) noexcept;

	int trainLSSVM(float gamma, int kernel,
		float * kernelParams, int kernelParamsCount,
		int modelParamsCount, int precomputedWeightsCount,
		float * dataset, int datasetRowCount, int datasetColCount,
		float * labels, int labelsLength,
		float * normalizedData, int normalizedDataRowCount, int normalizedDataColCount,
		float * means, int meansLength,
		float * stdDevs, int stdDevsLength,
		float * modelParams, int modelParamsLength,
		float * precomputedWeights, int precomputedWeightsLength) noexcept;

	int predictByLSSVM(int kernel,
		float* kernelParams, int kernelParamsCount,
		float* normalizedData, int normalizedDataRowCount, int normalizedDataColCount,
		float* labels, int labelsLength,
		float* means, int meansLength,
		float* stdDevs, int stdDevsLength,
		float* modelParams, int modelParamsLength,
		float* precomputedWeights, int precomputedWeightsLength,
		float* targetData, int targetDataRowCount, int targetDataColCount,
		float* prediction, int predictionLength) noexcept;

	int trainAndAnalyzeLSSVM(float gamma, int kernel,
		float * kernelParams, int kernelParamsCount,
		int modelParamsCount, int precomputedWeightsCount,
		int confusionMatrixElementsCount,
		float * dataset, int datasetRowCount, int datasetColCount,
		float * labels, int labelsLength,
		float * normalizedData, int normalizedDataRowCount, int normalizedDataColCount,
		float * means, int meansLength,
		float * stdDevs, int stdDevsLength,
		float * modelParams, int modelParamsLength,
		float * precomputedWeights, int precomputedWeightsLength,
		float * predictedLabels, int predictedLabelsLength,
		float* correctness, int correctnessLength,
		int* consfusionMatrix, int consfusionMatrixLength) noexcept;
}

//name: generateDataset
//input: int kernel
//input: column kernelParams
//input: int samplesCount
//input: int featuresCount
//input: double min 
//input: double max
//input: double violatorsPercentage
//output: column_list dataset [new(samplesCount, featuresCount)]
//output: column labels [new(samplesCount)]
EMSCRIPTEN_KEEPALIVE
int generateDataset(int kernel,
	float * kernelParams, int kernelParamsCount,
	int samplesCount, int featuresCount,
	float min, float max,
	float violatorsPercentage,
	float * dataset, int rowCount, int colCount,
	float * labels, int labelsLength) noexcept	
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

//name: normalizeDataset
//input: column_list data
//output: column_list normalizedData [new(data.columnCount, data.rowCount)]
//output: column means [new(data.columnCount)]
//output: column stdDevs [new(data.columnCount)]
EMSCRIPTEN_KEEPALIVE
int normalizeDataset(float * dataset, int datasetRowCount, int datasetColCount,
	float * normalizedData, int normalizedDataRowCount, int normalizedDataColCount,
	float * means, int meansLength,
	float * stdDevs, int stdDevsLength)	noexcept
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

//name: trainLSSVM
//input: double gamma
//input: int kernel
//input: column kernelParams
//input: int modelParamsCount
//input: int precomputedWeightsCount
//input: column_list dataset
//input: column labels
//output: column_list normalizedData [new(dataset.columnCount, dataset.rowCount)]
//output: column means [new(dataset.columnCount)]
//output: column stdDevs [new(dataset.columnCount)]
//output: column modelParams [new(modelParamsCount)]
//output: column precomputedWeights [new(precomputedWeightsCount)]
EMSCRIPTEN_KEEPALIVE
int trainLSSVM(float gamma, int kernel,
	float * kernelParams, int kernelParamsCount,
	int modelParamsCount, int precomputedWeightsCount,
	float * dataset, int datasetRowCount, int datasetColCount,
	float * labels, int labelsLength,
	float * normalizedData, int normalizedDataRowCount, int normalizedDataColCount,
	float * means, int meansLength,
	float * stdDevs, int stdDevsLength,
	float * modelParams, int modelParamsLength,
	float * precomputedWeights, int precomputedWeightsLength) noexcept
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

//name: predictByLSSVM
//input: int kernel
//input: column kernelParams
//input: column_list normalizedData
//input: column labels
//input: column means
//input: column stdDevs
//input: column modelParams
//input: column precomputedWeights
//input: column_list targetData
//output: column prediction [new(targetData.rowCount)]
EMSCRIPTEN_KEEPALIVE
int predictByLSSVM(int kernel,
	float * kernelParams, int kernelParamsCount,
	float * normalizedData, int normalizedDataRowCount, int normalizedDataColCount,
	float * labels, int labelsLength,
	float * means, int meansLength,
	float * stdDevs, int stdDevsLength,
	float * modelParams, int modelParamsLength,
	float * precomputedWeights, int precomputedWeightsLength,
	float * targetData, int targetDataRowCount, int targetDataColCount,
	float * prediction, int predictionLength) noexcept
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

//name: trainAndAnalyzeLSSVM
//input: double gamma
//input: int kernel
//input: column kernelParams
//input: int modelParamsCount
//input: int precomputedWeightsCount
//input: int confusionMatrixElementsCount
//input: column_list dataset
//input: column labels
//output: column_list normalizedData [new(dataset.columnCount, dataset.rowCount)]
//output: column means [new(dataset.columnCount)]
//output: column stdDevs [new(dataset.columnCount)]
//output: column modelParams [new(modelParamsCount)]
//output: column precomputedWeights [new(precomputedWeightsCount)]
//output: column predictedLabels [new(dataset.rowCount)]
//output: column correctness [new(dataset.rowCount)]
//output: column consfusionMatrix [new(confusionMatrixElementsCount)]
EMSCRIPTEN_KEEPALIVE
int trainAndAnalyzeLSSVM(float gamma, int kernel,
	float * kernelParams, int kernelParamsCount,
	int modelParamsCount, int precomputedWeightsCount,
	int confusionMatrixElementsCount,
	float * dataset, int datasetRowCount, int datasetColCount,
	float * labels, int labelsLength,
	float * normalizedData, int normalizedDataRowCount, int normalizedDataColCount,
	float * means, int meansLength,
	float * stdDevs, int stdDevsLength,
	float * modelParams, int modelParamsLength,
	float * precomputedWeights, int precomputedWeightsLength,
	float * predictedLabels, int predictedLabelsLength,
	float * correctness, int correctnessLength,
	int * consfusionMatrix, int consfusionMatrixLength) noexcept
{
	using namespace svm;

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
		|| (precomputedWeightsLength != datasetColCount + 1)
		|| (predictedLabelsLength != labelsLength)
		|| (correctnessLength != labelsLength)
		|| (consfusionMatrixLength != dmt::CONFUSION_MATR_SIZE))
		return INCORRECT_SIZE;

	// normalize data
	int resCode = dmt::getNormalizedDataset(dataset, datasetRowCount, datasetColCount,
		normalizedData, means, stdDevs);
	if (resCode != dmt::NO_ERRORS)
		return resCode;

	// train LS-SVM model
	resCode = trainLSSVM(gamma, kernel, kernelParams,
		normalizedData, labels, datasetRowCount, datasetColCount,
		modelParams, precomputedWeights);
	if (resCode != NO_ERRORS)
		return resCode;

	// get predictions
	resCode = predictByLSSVM(kernel, kernelParams, normalizedData, labels,
		normalizedDataColCount, normalizedDataRowCount,
		means, stdDevs, modelParams, precomputedWeights,
		dataset, predictedLabels, datasetRowCount);
	if (resCode != NO_ERRORS)
		return resCode;

	// analyze results
	return dmt::compareLabelsAndTheirPredictions(labels, predictedLabels, correctness,
		datasetRowCount, consfusionMatrix);
} // trainAndAnalyzeLSSVM
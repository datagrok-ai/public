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
		float * labels, int labelsLength);		

	int normalizeDataset(float * dataset, int datasetRowCount, int datasetColCount,
		float * normalizedData, int normalizedDataRowCount, int normalizedDataColCount,
		float * means, int meansLength,
		float * stdDevs, int stdDevsLength);
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
	float * labels, int labelsLength)	
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
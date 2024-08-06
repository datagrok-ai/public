#include <emscripten.h>

// The following provides convenient naming of the exported functions.
extern "C" {

	int train(float* features, int samplesCount, int featuresCount,
		float missingValue,
		float* labels, int labelsCount,
		int iterationsCount, float eta, int maxDepth, float lambda, float alpha,
		int* modelSize, int modelSizeCount,
		int* modelBytes, int modelBytesCount) noexcept;

	int predict(float* features, int samplesCount, int featuresCount,
		float missingValue,
		int* modelBytes, int modelBytesCount,
		float* predictions, int predictionsCount) noexcept;
}

#include "xgboost/include/xgboost/c_api.h"

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

EMSCRIPTEN_KEEPALIVE
int train(float * features, int samplesCount, int featuresCount,
	float missingValue,
	float * labels, int labelsCount,
	int iterationsCount, float eta, int maxDepth, float lambda, float alpha,
	int * modelSize, int modelSizeCount,
	int * modelBytes, int modelBytesCount) noexcept {
	if (labelsCount != samplesCount)
		return -1;

	// 0. Transform col-by-col features matrix to row-by-row
	Map<Matrix<float, Dynamic, Dynamic, ColMajor>> columnDataMatrix(features, samplesCount, featuresCount);
	Matrix<float, Dynamic, Dynamic, RowMajor> dataMatrix = columnDataMatrix;
	float* dataPtr = &dataMatrix(0);
	
	// 1. XGBoost data matrix
	DMatrixHandle dmatrix;
	XGDMatrixCreateFromMat(dataPtr, samplesCount, featuresCount, missingValue, &dmatrix);

	// 2. XGBooster labels
	XGDMatrixSetFloatInfo(dmatrix, "label", labels, samplesCount);

	// 3. Create model & specify settings
	DMatrixHandle eval_dmats[1] = { dmatrix };
	BoosterHandle booster;

	// settings as strings
	auto etaStr = std::to_string(eta);
	auto maxDepthStr = std::to_string(maxDepth);
	auto lambdaStr = std::to_string(lambda);
	auto alphaStr = std::to_string(alpha);	

	// specify settings
	XGBoosterCreate(eval_dmats, 1, &booster);
	XGBoosterSetParam(booster, "booster", "gbtree");
	XGBoosterSetParam(booster, "eta", etaStr.c_str());
	XGBoosterSetParam(booster, "max_depth", maxDepthStr.c_str());
	XGBoosterSetParam(booster, "lambda", lambdaStr.c_str());
	XGBoosterSetParam(booster, "alpha", alphaStr.c_str());

	// 4. Train model
	for (int i = 0; i < iterationsCount; ++i)
		XGBoosterUpdateOneIter(booster, i, dmatrix);

	// 5. Pack model
	char const packConfig[] = "{\"format\": \"json\"}";
	uint64_t outLen;
	char const* outDptr = NULL;
	XGBoosterSaveModelToBuffer(booster, packConfig, &outLen, &outDptr);

	// 6. Store model bytes

	// check reserved bytes
	if (sizeof(int) * modelBytesCount < outLen) {
		modelSize[0] = -1;
		return -1;
	}

	// store the required bytes count
	modelSize[0] = outLen;

	// copy bytes (it's recommended in the XGBoost lib)
	std::memcpy(reinterpret_cast<void*>(modelBytes), reinterpret_cast<const void*>(outDptr), outLen);

	
	// 7. Clearing
	XGDMatrixFree(dmatrix);
	XGBoosterFree(booster);

	return 0;
} // train

EMSCRIPTEN_KEEPALIVE
int predict(float* features, int samplesCount, int featuresCount,
	float missingValue,
	int* modelBytes, int modelBytesCount,
	float* predictions, int predictionsCount) noexcept {

	// 0. Transform col-by-col features matrix to row-by-row
	Map<Matrix<float, Dynamic, Dynamic, ColMajor>> columnDataMatrix(features, samplesCount, featuresCount);
	Matrix<float, Dynamic, Dynamic, RowMajor> dataMatrix = columnDataMatrix;
	float* dataPtr = &dataMatrix(0);

	// 1. XGBoost data matrix
	DMatrixHandle dmatrix;
	XGDMatrixCreateFromMat(dataPtr, samplesCount, featuresCount, missingValue, &dmatrix);

	// 2. Unpack model
	BoosterHandle booster;
	XGBoosterCreate(NULL, 0, &booster);
	XGBoosterLoadModelFromBuffer(booster, reinterpret_cast<void*>(modelBytes), modelBytesCount);

	// 3. Apply model
	uint64_t const* outShape;
	uint64_t outDim;
	char const config[] =
		"{\"training\": false, \"type\": 0, "
		"\"iteration_begin\": 0, \"iteration_end\": 0, \"strict_shape\": false}";	
	float const* outResult = NULL;
	XGBoosterPredictFromDMatrix(booster, dmatrix, config, &outShape, &outDim, &outResult);

	// 4. Store predictions
	for (int i = 0; i < samplesCount; ++i)
		predictions[i] = outResult[i];

	// 5. Clearing
	XGDMatrixFree(dmatrix);
	XGBoosterFree(booster);

	return 0;
} // predict

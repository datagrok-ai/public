#include "regression.h"

#include <emscripten.h>

// The following provides convenient naming of the exported functions.
extern "C" {
	int fitLinearRegressionParamsWithDataNormalizing(
	  float* features, int featuresRowCount, int featuresColCount,
	  float* featureAvgs, int featureAvgsCount,
	  float* featureStdDevs, int featureStdDevsCount,
	  float* targets, int targetsCount,
	  float targetsAvg,
	  float targetsStdDev,
	  float* params, int paramsCount) noexcept;

	int fitLinearRegressionParams(
	  float* features, int featuresRowCount, int featuresColCount,
	  float* targets, int targetsCount,
	  float* params, int paramsCount) noexcept;
}

//name: fitLinearRegressionParamsWithDataNormalizing
//input: column_list features
//input: column featureAvgs
//input: column featureStdDevs
//input: column targets
//input: double targetsAvg
//input: double targetsStdDev
//output: column params [new(features.rowCount + 1)]
EMSCRIPTEN_KEEPALIVE
int fitLinearRegressionParamsWithDataNormalizing(
	float* features, int featuresRowCount, int featuresColCount,
	float* featureAvgs, int featureAvgsCount,
	float* featureStdDevs, int featureStdDevsCount,
	float* targets, int targetsCount,
	float targetsAvg,
	float targetsStdDev,
	float* params, int paramsCount) noexcept {
	// Check sizes
	if ((featuresRowCount != targetsCount) || (featuresColCount != paramsCount - 1))
		return regr::INCORRECT_SIZES;

	return regr::fitLinearRegressionParams(features, featureAvgs, featureStdDevs, targets, targetsAvg, targetsStdDev, params, targetsCount, featuresColCount);
}

//name: fitLinearRegressionParams
//input: column_list features
//input: column targets
//output: column params [new(features.rowCount + 1)]
EMSCRIPTEN_KEEPALIVE
int fitLinearRegressionParams(
	float* features, int featuresRowCount, int featuresColCount,
	float* targets, int targetsCount,
	float* params, int paramsCount) noexcept {
	// Check sizes
	if ((featuresRowCount != targetsCount) || (featuresColCount != paramsCount - 1))
		return regr::INCORRECT_SIZES;

	return regr::fitLinearRegressionParams(features, targets, params, targetsCount, featuresColCount);
}
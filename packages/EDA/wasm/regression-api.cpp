#include "regression.h"

#include <emscripten.h>

// The following provides convenient naming of the exported functions.
extern "C" {
	int fitLinearRegressionParamsWithDataNormalizing(
	  float * features, int featuresRowCount, int featuresColCount,
	  float * featureAvgs, int featureAvgsCount,
	  float * featureStdDevs, int featureStdDevsCount,
	  float * targets, int targetsCount,
	  float targetsAvg,
	  float targetsStdDev,
	  int paramsCount,
	  float * params, int paramsLength) noexcept;

	int fitLinearRegressionParams(
	  float * features, int featuresRowCount, int featuresColCount,
	  float * targets, int targetsCount,
	  int paramsCount,
	  float * params, int paramsLength) noexcept;
}

//name: fitLinearRegressionParamsWithDataNormalizing
//input: column_list features
//input: column featureAvgs
//input: column featureStdDevs
//input: column targets
//input: double targetsAvg
//input: double targetsStdDev
//input: int paramsCount
//output: column params [new(paramsCount)]
EMSCRIPTEN_KEEPALIVE
int fitLinearRegressionParamsWithDataNormalizing(
	float * features, int featuresRowCount, int featuresColCount,
	float * featureAvgs, int featureAvgsCount,
	float * featureStdDevs, int featureStdDevsCount,
	float * targets, int targetsCount,
	float targetsAvg,
	float targetsStdDev,
	int paramsCount,
	float * params, int paramsLength) noexcept {
	// Check sizes
	if ((featuresRowCount != targetsCount) || (featuresColCount != paramsLength - 1))
		return regr::INCORRECT_SIZES;

	return regr::fitLinearRegressionParams(features, featureAvgs, featureStdDevs, targets, targetsAvg, targetsStdDev, params, targetsCount, featuresColCount);
}

//name: fitLinearRegressionParams
//input: column_list features
//input: column targets
//input: int paramsCount
//output: column params [new(paramsCount)]
EMSCRIPTEN_KEEPALIVE
int fitLinearRegressionParams(
	float * features, int featuresRowCount, int featuresColCount,
	float * targets, int targetsCount,
	int paramsCount,
	float * params, int paramsLength) noexcept {
	// Check sizes
	if ((featuresRowCount != targetsCount) || (featuresColCount != paramsLength - 1))
		return regr::INCORRECT_SIZES;

	return regr::fitLinearRegressionParams(features, targets, params, targetsCount, featuresColCount);
}
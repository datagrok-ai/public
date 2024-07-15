#ifndef XGBOOST_API_H
#define XGBOOST_API_H

int train(float* features, int samplesCount, int featuresCount,
	float missingValue,
	float* labels, int labelsCount,
	int iterationsCount, float eta, int maxDepth, float lambda, float alpha,
	int* modelSize, int modelSizeCount,
	int* modelBytes, int modelBytesCount);

int predict(float* features, int samplesCount, int featuresCount,
	float missingValue, 
	int* modelBytes, int modelBytesCount,
	int * predictions, int predictionsCount);

#endif // !XGBOOST_API_H


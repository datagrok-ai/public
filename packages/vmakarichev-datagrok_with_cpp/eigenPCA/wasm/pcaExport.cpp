#include <emscripten.h>

extern "C" {
    int principalComponentAnalysis(float * data,
	int height,
	int width,
	int numOfPrincipalComponents,
	float * principalComponents);

	int pcaWithApproximation(float * data,
	int height,
	int width,
	int numOfPrincipalComponents,
	float * principalComponents,
	float * approximation);

	float error(float * data1, float * data2, int length);
}

#include "PCA\PCA.h"

EMSCRIPTEN_KEEPALIVE
int principalComponentAnalysis(float * data,
	int height,
	int width,
	int numOfPrincipalComponents,
	float * principalComponents)
{
	return pca::pcaUsingCorrelationMatrix(data, height, width, numOfPrincipalComponents, principalComponents, 0);
}

EMSCRIPTEN_KEEPALIVE
int pcaWithApproximation(float * data,
	int height,
	int width,
	int numOfPrincipalComponents,
	float * principalComponents,
	float * approximation)
{
	return pca::pcaUsingCorrelationMatrix(data, height, width, numOfPrincipalComponents, principalComponents, approximation);
} 

EMSCRIPTEN_KEEPALIVE
float error(float * data1, float * data2, int length)
{
	return pca::mad(data1, data2, length);
} 



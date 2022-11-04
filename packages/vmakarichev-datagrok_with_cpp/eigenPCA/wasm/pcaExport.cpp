#include <emscripten.h>

extern "C" {
    int principalComponentAnalysis(float * data,
	int height,
	int width,
	int numOfPrincipalComponents,
	float * principalComponents,
	float * approxData = 0);
}

#include "PCA\PCA.h"
//#include "PCA\pca.cpp"

EMSCRIPTEN_KEEPALIVE
int principalComponentAnalysis(float * data,
	int height,
	int width,
	int numOfPrincipalComponents,
	float * principalComponents,
	float * approxData)
{
	return pca::pcaUsingCorrelationMatrix(data, height, width, numOfPrincipalComponents, principalComponents, approxData);
} // principalComponentAnalysis
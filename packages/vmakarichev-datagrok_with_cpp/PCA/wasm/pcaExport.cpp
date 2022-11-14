#include <emscripten.h>

extern "C" {
    int principalComponentAnalysis(float * data,
	      int dataNumOfRows,
	      int dataNumOfColumns,
	      int numOfPrincipalComponents,
	      float * principalComponents,
	      int principalComponentsNumOfRows,
	      int principalComponentsNumOfColumns);

	int pcaWithApproximation(float * data,
	      int dataNumOfRows,
	      int dataNumOfColumns,
	      int numOfPrincipalComponents,
	      float * principalComponents,
	      int principalComponentsNumOfRows,
	      int principalComponentsNumOfColumns,
	      float * approximation,
	      int approximationComponentsNumOfRows,
	      int approximationComponentsNumOfColumns);

	float error(float * data1, int data1Length, float * data2, int data2Length);
}

#include "PCA\PCA.h"

EMSCRIPTEN_KEEPALIVE
int principalComponentAnalysis(float * data,
      int dataNumOfRows,
	  int dataNumOfColumns,
	  int numOfPrincipalComponents,
	  float * principalComponents,
	  int principalComponentsNumOfRows,
	  int principalComponentsNumOfColumns)
{
	return pca::pcaUsingCorrelationMatrix(data, dataNumOfColumns, dataNumOfRows, 
	               numOfPrincipalComponents, principalComponents, 0);
}

EMSCRIPTEN_KEEPALIVE
int pcaWithApproximation(float * data,
	  int dataNumOfRows,
	  int dataNumOfColumns,
	  int numOfPrincipalComponents,
	  float * principalComponents,
	  int principalComponentsNumOfRows,
	  int principalComponentsNumOfColumns,
	  float * approximation,
	  int approximationComponentsNumOfRows,
	  int approximationComponentsNumOfColumns)
{
	return pca::pcaUsingCorrelationMatrix(data, dataNumOfColumns, dataNumOfRows, 
	               numOfPrincipalComponents, principalComponents, approximation);
} 

EMSCRIPTEN_KEEPALIVE
float error(float * data1, int data1Length, float * data2, int data2Length)
{
	return pca::mad(data1, data2, data1Length);
} 



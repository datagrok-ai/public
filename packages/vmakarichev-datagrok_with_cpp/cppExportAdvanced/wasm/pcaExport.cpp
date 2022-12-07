#include <emscripten.h>

extern "C" {
	
    int principalComponentAnalysis(float * data,
	      int dataNumOfRows,
	      int dataNumOfColumns,
	      int numOfPrincipalComponents,
	      float * principalComponents,
	      int principalComponentsNumOfRows,
	      int principalComponentsNumOfColumns);

	float error(float * data1, int data1Length, float * data2, int data2Length);
}

#include "PCA\PCA.h"

//top-menu: Tools | Data Science | PCA by Eigen
//tag: ML
//name: pca
//input: dataframe table
//input: column_list columns
//input: int componentsCount
//output: column_list components [new(columns.rowCount, componentsCount)]
//output: dataframe result [components]
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

//name: error
//input: dataframe df
//input: column col1
//input: column col2
//output: double mad [_callResult]
EMSCRIPTEN_KEEPALIVE
float error(float * data1, int data1Length, float * data2, int data2Length)
{
	return pca::mad(data1, data2, data1Length);
} 



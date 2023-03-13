// This file contains C++-functions that are exported to wasm.

// The tool Emscripten is applied (the header emscripten.h is included
// and each exported function is marked by EMSCRIPTEN_KEEPALIVE).

// Also, each function has a special DATAGROK annotation for C++-functions.
// This approach provides further usage of C++-to-wasm export script that 
// performes all routine steps.

#include <emscripten.h>

// The following provides convenient naming of the exported functions.
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

//top-menu: Tools | Data Science | PCA
//name: PCA
//input: dataframe table
//input: column_list columns
//input: int componentsCount
//output: column_list components [new(columns.rowCount, componentsCount)]
//output: dataframe PCA [components]
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



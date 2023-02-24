// This file contains C++-functions that are exported to wasm.

// The tool Emscripten is applied (the header emscripten.h is included
// and each exported function is marked by EMSCRIPTEN_KEEPALIVE).

// Also, each function has a special DATAGROK annotation for C++-functions.
// This approach provides further usage of C++-to-wasm export script that 
// performes all routine steps.

#include <emscripten.h>

// The following provides convenient naming of the exported functions.
extern "C" {

	int partialLeastSquareRegression(float * predictorColumns,
	       int rowCount,
		   int columnCount,
		   float * responseColumn,
		   int responceColumnLength,
		   int componentsCount,
		   float * predictionColumn,
		   int predictionColumnLength,
		   float * regressionCoefficients,
		   int regressionCoefficientsLength,
		   float * predictorScoresColumns,
		   int predictorScoresColumnsRowCount,
		   int predictorScoresColumnsColumnCount,
		   float * predictionScoresColumns,
		   int predictionScoresColumnsRowCount,
		   int predictionScoresColumnsColumnCount,
	       float * predictionLoadingsColumns,
	       int predictionLoadingsColumnsRowCount,
	       int predictionLoadingsColumnsColumnCount);    
}

#include "PLS\PLS.h"

//name: pls
//input: dataframe table
//input: column_list features
//input: column predict
//input: int componentsCount
//output: column prediction [new(predict.rowCount)]
//output: column regressionCoefficients [new(features.columnCount)]
//output: column_list tScores [new(predict.rowCount, componentsCount)]
//output: column_list uScores [new(predict.rowCount, componentsCount)]
//output: column_list xLoadings [new(features.columnCount, componentsCount)]
EMSCRIPTEN_KEEPALIVE
int partialLeastSquareRegression(float * predictorColumns,
	   int rowCount,
	   int columnCount,
	   float * responseColumn,
	   int responceColumnLength,
	   int componentsCount,
	   float * predictionColumn,
	   int predictionColumnLength,
	   float * regressionCoefficients,
	   int regressionCoefficientsLength,
	   float * predictorScoresColumns,
	   int predictorScoresColumnsRowCount,
	   int predictorScoresColumnsColumnCount,
	   float * predictionScoresColumns,
	   int predictionScoresColumnsRowCount,
	   int predictionScoresColumnsColumnCount,
	   float * predictionLoadingsColumns,
	   int predictionLoadingsColumnsRowCount,
	   int predictionLoadingsColumnsColumnCount)
{
	return pls::partialLeastSquareExtended(predictorColumns, rowCount, columnCount,
		responseColumn, componentsCount, predictionColumn, regressionCoefficients,
		predictorScoresColumns, predictionScoresColumns, predictionLoadingsColumns);
} 




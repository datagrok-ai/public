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
	       int predictionLoadingsColumnsColumnCount,
		   float * yLoadingsColumn,
		   int yLoadingsColumnLength);    
}

#include "PLS\PLS.h"

//name: partialLeastSquareRegression
//input: dataframe table
//input: column_list features
//input: column predict
//input: int componentsCount
//output: column prediction [new(predict.rowCount)]
//output: column regressionCoefficients [new(features.columnCount)]
//output: column_list tScores [new(predict.rowCount, componentsCount)]
//output: column_list uScores [new(predict.rowCount, componentsCount)]
//output: column_list xLoadings [new(features.columnCount, componentsCount)]
//output: column yLoadings [new(componentsCount)]
EMSCRIPTEN_KEEPALIVE
int partialLeastSquareRegression(float * featuresColumns,
	   int rowCount,
	   int columnCount,
	   float * predictColumn,
	   int predictColumnLength,
	   int componentsCount,
	   float * predictionColumn,
	   int predictionColumnLength,
	   float * regressionCoefficients,
	   int regressionCoefficientsLength,
	   float * tScoresColumns,
	   int tScoresColumnsRowCount,
	   int tScoresColumnsColumnCount,
	   float * uScoresColumns,
	   int uScoresColumnsRowCount,
	   int uScoresColumnsColumnCount,
	   float * xLoadingsColumns,
	   int xLoadingsColumnsRowCount,
	   int xLoadingsColumnsColumnCount,
	   float * yLoadingsColumn,
	   int yLoadingsColumnLength)
{
	return pls::partialLeastSquareExtended(featuresColumns, rowCount, columnCount,
		predictColumn, componentsCount, predictionColumn, regressionCoefficients,
		tScoresColumns, uScoresColumns, xLoadingsColumns, yLoadingsColumn);
} 




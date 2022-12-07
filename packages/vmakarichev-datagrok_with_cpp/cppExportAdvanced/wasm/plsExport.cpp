#include <emscripten.h>

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

//top-menu: Tools | Data Science | MVA (PLS) by Eigen
//name: pls
//tag: ML
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




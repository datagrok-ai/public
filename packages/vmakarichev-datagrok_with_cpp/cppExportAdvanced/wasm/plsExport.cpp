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
//input: column prediction [new(predict.rowCount)]
//input: column regressionCoefficients [new(features.columnCount)]
//input: column_list tScores [new(predict.rowCount, componentsCount)]
//input: column_list uScores [new(predict.rowCount, componentsCount)]
//input: column_list xLoadings [new(features.columnCount, componentsCount)]
//output: objects obj [prediction, regressionCoefficients, tScores, uScores, xLoadings]
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




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
		   int regressionCoefficientsLength);

	int partialLeastSquareRegressionExtended(float * predictorColumns,
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
		   int predictionScoresColumnsColumnCount);    
}

#include "PLS\PLS.h"

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
	   int regressionCoefficientsLength)
{
	return pls::partialLeastSquare(predictorColumns, rowCount, columnCount,
	               responseColumn, componentsCount, predictionColumn, regressionCoefficients);
}

EMSCRIPTEN_KEEPALIVE
int partialLeastSquareRegressionExtended(float * predictorColumns,
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
	   int predictionScoresColumnsColumnCount)
{
	return pls::partialLeastSquareExtended(predictorColumns, rowCount, columnCount,
		responseColumn, componentsCount, predictionColumn, regressionCoefficients,
		predictorScoresColumns, predictionScoresColumns);
} 




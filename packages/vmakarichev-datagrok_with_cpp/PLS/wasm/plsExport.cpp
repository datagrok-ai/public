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




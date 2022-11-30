// PLS.h
// Declarations of functions that provides Partial Least Square (PLS) Regression.

// An implementation of the algorithm PLS1 without X-deflation is used.
// Source paper: Ulf G. Indahl, The geometry of PLS1 explained properly: 
//               10 key notes on mathematical properties of and some alternative 
//               algorithmic approaches to PLS1 modelling, DOI: https://doi.org/10.1002/cem.2589
// Also, the following aricle is used: https://doi.org/10.1016/S0169-7439(01)00155-1

#ifndef PLS_H
#define PLS_H

namespace pls {

	typedef float Float;	
	typedef double Double;

	enum ResultCode { NO_ERROR = 0, UNCORRECT_ARGUMENTS_ERROR, COMPUTATION_ERROR, METHOD_ERROR };

	/* Partial Least Square (PLS1).
	      predictorColumnsDataPtr - data from columns that are used for prediction
		  rowCount - number of rows
		  columnCount - number of columns
		  responseColumnDataPtr - data from column that is predicted, i.e. responce
		  componentsCount - number of components that extracted in PLS
		  predictionDataPtr - prediction obtained using PLS (its size is equal to the size of responce)
		  regressionCoefficients - coeffcient of linear regression that are computed (their size is eqaul to the number of columns)
	*/
	int partialLeastSquare(Float * predictorColumnsDataPtr,
		const int rowCount,
		const int columnCount,
		Float * responseColumnDataPtr,
		const int componentsCount,
		Float * predictionDataPtr,
		Float * regressionCoefficients) noexcept;	


	/* Partial Least Square (PLS1) - extended version: scores data is provided.
	      predictorColumnsDataPtr - data from columns that are used for prediction (X)
	      rowCount - number of rows
	      columnCount - number of columns
	      responseColumnDataPtr - data from column that is predicted, i.e. responce (Y)
	      componentsCount - number of components that extracted in PLS (A)
	      predictionDataPtr - prediction obtained using PLS (its size is equal to the size of responce)
	      regressionCoefficientsPtr - coeffcient of linear regression that are computed (their size is eqaul to the number of columns) (b)
		  predictorScoresPtr - scores of predectors (T)
		  responceScoresPtr - scores of response (U)
		  predictorLoadingsPtr - loadings of predictors (P)
	*/
	int partialLeastSquareExtended(Float * predictorColumnsDataPtr,
		const int rowCount,
		const int columnCount,
		Float * responseColumnDataPtr,
		const int componentsCount,
		Float * predictionDataPtr,
		Float * regressionCoefficientsPtr,
		Float * predictorScoresPtr,
		Float * responceScoresPtr,
	    Float * predictorLoadingsPtr) noexcept;	
};

#endif



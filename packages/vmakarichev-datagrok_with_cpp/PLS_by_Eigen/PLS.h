// PLS.h
// Declarations of functions that provides Partial Least Square (PLS) Regression.

// An implementation of the algorithm PLS1 without X-deflation is used.
// Source paper: Ulf G. Indahl, The geometry of PLS1 explained properly: 
//               10 key notes on mathematical properties of and some alternative 
//               algorithmic approaches to PLS1 modelling, DOI: https://doi.org/10.1002/cem.2589

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
		  regressionCoefficients - coeffcient of linear regression that are computed (their size is eqaul to the number of columns).

	   WARNING! This provides correct results for the case, when predictor columns have the same distribution!	  
	*/
	int partialLeastSquare(Float * predictorColumnsDataPtr,
		const int rowCount,
		const int columnCount,
		Float * responseColumnDataPtr,
		const int componentsCount,
		Float * predictionDataPtr,
		Float * regressionCoefficients) noexcept;

	/* Partial Least Square (PLS1).
	      predictorColumnsDataPtr - data from columns that are used for prediction
	      rowCount - number of rows
	      columnCount - number of columns
	      responseColumnDataPtr - data from column that is predicted, i.e. responce
	      componentsCount - number of components that extracted in PLS
	      predictionDataPtr - prediction obtained using PLS (its size is equal to the size of responce)
	      regressionCoefficients - coeffcient of linear regression that are computed (their size is eqaul to the number of columns)

	WARNING! This provides correct results for the case, when predictor columns have the same distribution!

	This implementation of PLS uses deflation step, so it's slower!	*/
	int partialLeastSquare_slow(Float * predictorColumnsDataPtr,
		const int rowCount,
		const int columnCount,
		Float * responseColumnDataPtr,
		const int componentsCount,
		Float * predictionDataPtr,
		Float * regressionCoefficients) noexcept;

	/* Partial Least Square (PLS1).
	      predictorColumnsDataPtr - data from columns that are used for prediction
	      rowCount - number of rows
	      columnCount - number of columns
	      responseColumnDataPtr - data from column that is predicted, i.e. responce
	      componentsCount - number of components that extracted in PLS
	      predictionDataPtr - prediction obtained using PLS (its size is equal to the size of responce)
	      regressionCoefficients - coeffcient of linear regression that are computed (their size is eqaul to the number of columns)
	*/
	int partialLeastSquare_norm(Float * predictorColumnsDataPtr,
		const int rowCount,
		const int columnCount,
		Float * responseColumnDataPtr,
		const int componentsCount,
		Float * predictionDataPtr,
		Float * regressionCoefficients) noexcept;

	// Maximum absolute deviation between arrays
	Float mad(Float * arr1, Float * arr2, const int length) noexcept;
};

#endif



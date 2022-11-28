// pls.cpp
// Principal Component Analysis (PCA) using the lib Eigen: implementation of functions

#include<iostream>
using namespace std;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "PLS.h"
using pls::Float;
using pls::Double;

  /* Partial Least Square (PLS1).
  predictorColumnsDataPtr - data from columns that are used for prediction
  rowCount - number of rows
  columnCount - number of columns
  responseColumnDataPtr - data from column that is predicted, i.e. responce
  componentsCount - number of components that extracted in PLS
  predictionDataPtr - prediction obtained using PLS (its size is equal to the size of responce)
  regressionCoefficients - coeffcient of linear regression that are computed (their size is eqaul to the number of columns)
  */
int pls::partialLeastSquare(Float * predictorColumnsDataPtr,
	const int rowCount,
	const int columnCount,
	Float * responseColumnDataPtr,
	const int componentsCount,
	Float * predictionDataPtr,
	Float * regressionCoefficients) noexcept
{
	// check correctness of arguments
	if (componentsCount <= 0 || componentsCount > columnCount)
		return UNCORRECT_ARGUMENTS_ERROR;

	// Further, notation from the paper https://doi.org/10.1002/cem.2589 is used (see Algorithm 2).

	// create matrix, which is associated with predictor data
	Map < Matrix<Float, Dynamic, Dynamic, ColMajor>> D(predictorColumnsDataPtr, rowCount, columnCount);

	// compute mean value of each column of D
	Vector<Float, Dynamic> mu = D.colwise().mean();

	// mean-centered version of D
	Matrix<Float, Dynamic, Dynamic, ColMajor> X = D.rowwise() - mu.transpose();

    // vector for standard deviations of X
	Vector<Float, Dynamic> stdDevX(columnCount);
	
	Float rowCountSqrt = sqrt(static_cast<Float>(rowCount));

    // normilizing X-columns
	for (int i = 0; i < columnCount; i++)
	{
		stdDevX(i) = X.col(i).norm() / rowCountSqrt;
		X.col(i) = X.col(i) / stdDevX(i);		
	}	

	// create a vector, which is associated with responce or predicted data 
	Map<Vector<Float, Dynamic>> ySource(responseColumnDataPtr, rowCount);		

    // mean value of the responce
	Vector<Float, 1> meanY;
	meanY(0) = ySource.mean();		

    // mean-centered version of the responce
	Vector<Float, Dynamic> y = ySource.rowwise() - meanY;

    // standard deviation
	Float stdDevY = sqrt(y.squaredNorm() / rowCount);

    // normalizing
	y /= stdDevY;	

	// create a vector, which is associtated with regression coefficients
	Map<Vector<Float, Dynamic>> b(regressionCoefficients, columnCount);

	// create a vector, which is associated with prediction data 
	Map<Vector<Float, Dynamic>> prediction(predictionDataPtr, rowCount);

	// PLS1 algorithm routine

	Matrix<Float, Dynamic, Dynamic, ColMajor> W(columnCount, componentsCount);

	Matrix<Float, Dynamic, Dynamic, ColMajor> P(columnCount, componentsCount);

	Matrix<Float, Dynamic, Dynamic, ColMajor> T(rowCount, componentsCount);

	Vector<Float, Dynamic> normTau(componentsCount);

	Vector<Float, Dynamic> q(componentsCount);

	Vector<Float, Dynamic> normV(componentsCount);

	// PLS1 algorithm: see Algorithm 2 in https://doi.org/10.1002/cem.2589

	Vector<Float, Dynamic> w = (X.transpose() * y);

	normV(0) = w.norm();

	// prevent division by zero
	if (normV(0) == static_cast<Float>(0))
		return METHOD_ERROR;

	w = w / normV(0);

	W.col(0) = w;

	Vector<Float, Dynamic> t = X * w;

	normTau(0) = t.norm();

	// prevent division by zero
	if (normTau(0) == static_cast<Float>(0))
		return METHOD_ERROR;

	t = t / normTau(0);

	T.col(0) = t;

	Vector<Float, Dynamic> p = X.transpose() * t;

	P.col(0) = p;

	q(0) = t.transpose() * y;

	for (int a = 1; a < componentsCount; a++)
	{
		w = normV(a - 1) * (w - p / normTau(a - 1));

		normV(a) = w.norm();

		// prevent division by zero
		if (normV(a) == static_cast<Float>(0))
			return METHOD_ERROR;

		w = w / normV(a);

		W.col(a) = w;

		t = X * w;

		t = t - T.leftCols(a) * (T.leftCols(a).transpose() * t);

		normTau(a) = t.norm();

		// prevent division by zero
		if (normTau(a) == static_cast<Float>(0))
			return METHOD_ERROR;

		t = t / normTau(a);

		T.col(a) = t;

		p = X.transpose() * t;

		P.col(a) = p;

		q(a) = t.transpose() * y;
	} // for a	

	  // compute coefficients of regression
	Matrix<Float, Dynamic, Dynamic> H = P.transpose() * W;

	// chech existence of inverse matrix
	if (H.determinant() == static_cast<Float>(0))
		return METHOD_ERROR;

	b = W * H.inverse() * q;

	for (int i = 0; i < columnCount; i++)
		b(i) *= stdDevY / stdDevX(i);   
	
	prediction = D * b;	

	return NO_ERROR;
} // partialLeastSquare

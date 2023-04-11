// pls.cpp
// Principal Component Analysis (PCA) using the lib Eigen: implementation of functions

// The following STL lib is used for printing results and their verifying
//#include<iostream>
//using namespace std;

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
	
	// TODO: to discuss a constant term of the regression
	// a constant term
	//Vector<Float, 1> shift;
	//shift(0) = ySource(0) - D.row(0) * b;
	//q(0) - P.col(0).transpose().dot(b);	
	//prediction = (D * b).rowwise() + shift;

	prediction = D * b;	

	return NO_ERROR;
} // partialLeastSquare

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
int pls::partialLeastSquareExtended(Float * predictorColumnsDataPtr,
	const int rowCount,
	const int columnCount,
	Float * responseColumnDataPtr,
	const int componentsCount,
	Float * predictionDataPtr,
	Float * regressionCoefficientsPtr,
	Float * predictorScoresPtr,
	Float * responceScoresPtr,
	Float * predictorLoadingsPtr) noexcept
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
	
	// standard deviations of X
	Vector<Float, Dynamic> stdDevX(columnCount);
		
	Float rowCountSqrt = sqrt(static_cast<Float>(rowCount));

	// normalizing X
	for (int i = 0; i < columnCount; i++)
	{
		stdDevX(i) = X.col(i).norm() / rowCountSqrt;

        // check deviation
		if(stdDevX(i) == static_cast<Float>(0))
		    return UNCORRECT_ARGUMENTS_ERROR;

		X.col(i) = X.col(i) / stdDevX(i);
	}	

	// create a vector, which is associated with responce or predicted data 
	Map<Vector<Float, Dynamic>> ySource(responseColumnDataPtr, rowCount);

	// mean value of Y: Eigen vector is used in order to provide broadcasting
	Vector<Float, 1> meanY;
	meanY(0) = ySource.mean();

	// centering Y
	Vector<Float, Dynamic> y = ySource.rowwise() - meanY;

	// standard deviation of Y normalizing Y
	Float stdDevY = sqrt(y.squaredNorm() / rowCount);

	// check deviation
	if(stdDevY == static_cast<Float>(0))
	    return UNCORRECT_ARGUMENTS_ERROR;

	// normalizing Y
	y /= stdDevY;
	
	// create a vector, which is associtated with regression coefficients
	Map<Vector<Float, Dynamic>> b(regressionCoefficientsPtr, columnCount);

	// create a vector, which is associated with prediction data 
	Map<Vector<Float, Dynamic>> prediction(predictionDataPtr, rowCount);

	// weights matrix, W
	Matrix<Float, Dynamic, Dynamic, ColMajor> W(columnCount, componentsCount);

	// X-loadings matrix, P
	Map<Matrix<Float, Dynamic, Dynamic, ColMajor>> P(predictorLoadingsPtr, columnCount, componentsCount);

	//Matrix<Float, Dynamic, Dynamic, ColMajor> P(columnCount, componentsCount);

	// X-scores, T
	Map<Matrix<Float, Dynamic, Dynamic, ColMajor>> T(predictorScoresPtr, rowCount, componentsCount);

	// Y-scores, U
	Map<Matrix<Float, Dynamic, Dynamic, ColMajor>> U(responceScoresPtr, rowCount, componentsCount);

	// Y-loadings, q
	Vector<Float, Dynamic> q(componentsCount);

	// PLS1 routine auxiliry vectors
	Vector<Float, Dynamic> normTau(componentsCount);
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

	// compute Y-scores
	U = y * q.transpose() / q.squaredNorm();

	// compute coefficients of regression
	Matrix<Float, Dynamic, Dynamic> H = P.transpose() * W;

	// chech existence of inverse matrix
	if (H.determinant() == static_cast<Float>(0))
		return METHOD_ERROR;

	// compute regression coefficients
	b = W * H.inverse() * q;

	// ... also, we take into account a normalizing
	for (int i = 0; i < columnCount; i++)
		b(i) *= stdDevY / stdDevX(i);

	// compute predictions
	prediction = D * b;	

	// Remove the following comments in order to print and verify results
	//cout << "\nW_star:\n" << Wstar << endl;	
	//cout << "\nU:\n" << U << endl;
	//cout << "\nU.tr * U:\n" << U.transpose() * U << endl; // this must be identity matrix	
	//cout << "\nb:\n" << b << endl;
	//cout << "\nq:\n" << q << endl;
	//cout << "\nD:\n" << D << endl;	
	//cout << "\nP:\n" << P << endl;
	//cout << "\nT:\n" << T << endl;
	//cout << "\nT.tr * T:\n" << T.transpose() * T << endl; // this must be identity matrix
	//cout << "\nW:\n" << W << endl;
	//cout << "\nW.tr * W:\n" << W.transpose() * W << endl; // this must be identity matrix
	//cout << "\nprediction\n" << prediction << endl;

	return NO_ERROR;
} // partialLeastSquareExtended




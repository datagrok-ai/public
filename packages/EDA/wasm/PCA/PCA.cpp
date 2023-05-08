// PCA.cpp
// Principal Component Analysis using the lib Eigen: implementations of functions

#include "../../../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "PCA.h"
using pca::Float;
using pca::Integer;
using pca::Double;

/* Principal Component Analysis of the data: using correlation matrix.
     data - input matrix;
     height, width - sizes of the input;
     numOfPrincipalComponents - number of principal components to be computed;
     principalComponents - the principal components computed;
	 approxData - approximation of the input data using principal components obtained.  */
int pca::pcaUsingCorrelationMatrix(Float * data,
	const int height,
	const int width,
	const int numOfPrincipalComponents,
	const int centerNum,
	const int scaleNum,
	Float * principalComponents,
	Float * approxData) noexcept
{
	/* Here, we use a MODIFICATION of the algorithm given in 
	   Charu C. Aggarwal. Data Mining: The Textbook. Springer, 2015,
	   (see page 42). */

	// check number of principal components
	if (height < numOfPrincipalComponents || numOfPrincipalComponents < 1)
		return UNCORRECT_ARGUMENTS_ERROR;

	// assign data and Eigen matrix
	Map< Matrix<Float, Dynamic, Dynamic, RowMajor> > dataMatrix(data, height, width);

	Vector<Float, Dynamic> means = dataMatrix.rowwise().mean();		

	if (centerNum != 0)		
		dataMatrix = dataMatrix.colwise() - means;

	if (scaleNum != 0)
	    dataMatrix = dataMatrix.rowwise().normalized() * sqrt(height);
	
	Matrix<Float, Dynamic, Dynamic>	corMatrix = dataMatrix * dataMatrix.transpose();	

	// The following solver computes eigen vals & vectors: the order of eigen vals is increasing.
	SelfAdjointEigenSolver<Matrix<Float, Dynamic, Dynamic>> eigensolver(corMatrix);

	// Check result of eigen values & vectors computation.
	if (eigensolver.info() != Success)
		return COMPUTATION_ERROR;

	// Check order of computed eigen values: increasing order is expected
	Vector<Float, Dynamic> eigenVals = eigensolver.eigenvalues();
	for(int i = 1; i < eigenVals.size(); i++)
	    if(eigenVals(i - 1) > eigenVals(i))
		    return METHOD_ERROR;
	
	// get feature vectors, taking into account increasing order of computed eigen values
	Matrix<Float, Dynamic, Dynamic, ColMajor> featureVectors
		= (eigensolver.eigenvectors().rowwise().reverse())(all, seq(0, numOfPrincipalComponents - 1));

	// assign principal components and Eigen matrix
	Map< Matrix<Float, Dynamic, Dynamic, RowMajor> > 
	     princCompMatrix(principalComponents, numOfPrincipalComponents, width);

	princCompMatrix = featureVectors.transpose() * dataMatrix; 

	// computation of approximation
	if (approxData != NULL)
	{
		// assign data and Eigen matrix
		Map< Matrix<Float, Dynamic, Dynamic, RowMajor> > approxMatrix(approxData, height, width);

		approxMatrix = (featureVectors * princCompMatrix).colwise() + means;
	}	

	return NO_ERROR;
} // pcaUsingCorrelationMatrix

/*{
	// Here, we use a MODIFICATION of the algorithm given in 
	// Charu C. Aggarwal. Data Mining: The Textbook. Springer, 2015,
	// (see page 42). 

	// check number of principal components
	if (height < numOfPrincipalComponents || numOfPrincipalComponents < 1)
		return UNCORRECT_ARGUMENTS_ERROR;

	// assign data and Eigen matrix
	Map< Matrix<Float, Dynamic, Dynamic, RowMajor> > dataMatrix(data, height, width);
	
	Vector<Float, Dynamic> means = dataMatrix.rowwise().mean();

	Matrix<Float, Dynamic, Dynamic> corMatrix 
		= dataMatrix * dataMatrix.transpose() / width - means * means.transpose();

	// The following solver computes eigen vals & vectors: the order of eigen vals is increasing.
	SelfAdjointEigenSolver<Matrix<Float, Dynamic, Dynamic>> eigensolver(corMatrix);

	// Check result of eigen values & vectors computation.
	if (eigensolver.info() != Success)
		return COMPUTATION_ERROR;

	// Check order of computed eigen values: increasing order is expected
	Vector<Float, Dynamic> eigenVals = eigensolver.eigenvalues();
	for(int i = 1; i < eigenVals.size(); i++)
	    if(eigenVals(i - 1) > eigenVals(i))
		    return METHOD_ERROR;
	
	// get feature vectors, taking into account increasing order of computed eigen values
	Matrix<Float, Dynamic, Dynamic, ColMajor> featureVectors
		= (eigensolver.eigenvectors().rowwise().reverse())(all, seq(0, numOfPrincipalComponents - 1));

	// assign principal components and Eigen matrix
	Map< Matrix<Float, Dynamic, Dynamic, RowMajor> > 
	     princCompMatrix(principalComponents, numOfPrincipalComponents, width);

	// compute principal componets
	princCompMatrix = featureVectors.transpose() * (dataMatrix.colwise() - means);

	// computation of approximation
	if (approxData != NULL)
	{
		// assign data and Eigen matrix
		Map< Matrix<Float, Dynamic, Dynamic, RowMajor> > approxMatrix(approxData, height, width);

		approxMatrix = (featureVectors * princCompMatrix).colwise() + means;
	}	

	return NO_ERROR;
} */

// Maximum absolute deviation between arrays
Float pca::mad(Float * arr1, Float * arr2, const int length) noexcept
{
	// Solution using Eigen: nice, but additional structures are created! 
	/*Map<Vector<Float, Dynamic>> vec1(arr1, length);
	Map<Vector<Float, Dynamic>> vec2(arr2, length);
	return ((vec1 - vec2).cwiseAbs()).maxCoeff();*/

	// Naive solution
	Float result = fabs(arr1[0] - arr2[0]);

	for (int i = 1; i < length; i++)
		result = fmax(result, fabs(arr1[i] - arr2[i]));

	return result;
}
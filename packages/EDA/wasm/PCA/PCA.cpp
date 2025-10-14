// PCA.cpp
// Principal Component Analysis using the lib Eigen: implementations of functions

#include "../../../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "PCA.h"
using pca::Float;
using pca::Integer;
using pca::Double;
using pca::MAX_ITER;
using pca::TOL;

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

// Maximum absolute deviation between arrays
Float pca::mad(Float * arr1, Float * arr2, const int length) noexcept
{
	// Naive solution
	Float result = fabs(arr1[0] - arr2[0]);

	for (int i = 1; i < length; i++)
		result = fmax(result, fabs(arr1[i] - arr2[i]));

	return result;
}

/* The NIPALS algorithm for PCA.
     data - input matrix;
     height, width - sizes of the input;
	 numOfPrincipalComponents - number of principal components to be computed;
	 principalComponents - the principal components computed.

	 Reference
	 H. Risvik. Principal component analysis (PCA) & NIPALS algorithm.
	 https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=a1a75889929b604cfefcb7cddfd58accce8d01d8
*/
int pca::nipals(Float* data,
	const int height,
	const int width,
	const int numOfPrincipalComponents,
	Float* principalComponents) noexcept
{
	Map< Matrix<Float, Dynamic, Dynamic, ColMajor> > X(data, height, width);	
	Vector<Float, Dynamic> means = X.colwise().mean();
    X = X.rowwise() - means.transpose();
    X = X.colwise().normalized() * sqrt(height);

	Map<Matrix<Float, Dynamic, Dynamic, ColMajor>> T(principalComponents, height, numOfPrincipalComponents);	
	Vector<Float, Dynamic> p(width);
	Vector<Float, Dynamic> t(height);
	Vector<Float, Dynamic> tNew(height);

	for (int i = 0; i < numOfPrincipalComponents; ++i) {
		t = X.col(i);

		for (int j = 0; j < MAX_ITER; ++j) {
			p = X.transpose() * t / t.squaredNorm();
			p.normalize();
			tNew = X * p;

			if ((tNew - t).norm() < TOL)
				break;

			t = tNew;
		}

		T.col(i) = t;
		X = X - t * p.transpose();		
	}

	return 0;
}
// PCA.cpp
// Principal Component Analysis using the lib Eigen: implementations of functions

#include "..\..\..\Eigen\Eigen\Dense"
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
	Float * principalComponents,
	Float * approxData) noexcept
{
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
}

// Cumpute principal components of the data.
int computePrincipalComponents(void ** data,
	const Float * meanValsOfRows,
	const Float * featureVectors,
	const int heightOfFloats,
	const int heightOfInts,
	const int width,
	const int numOfPrincipalComponents,
	Float ** principalComponents) noexcept
{
	int dim = heightOfFloats + heightOfInts;

	// multiply float data
	for (int j = 0; j < numOfPrincipalComponents; j++)
		for (int i = 0; i < width; i++)
			principalComponents[j][i] = 0;

	for (int s = 0; s < heightOfFloats; s++)
	{
		Float * left = static_cast<Float *>(data[s]);
		Float mean = meanValsOfRows[s];

		for (int j = 0; j < numOfPrincipalComponents; j++)
		{
			Float r = featureVectors[s * numOfPrincipalComponents + j]; // (s, j);

			for (int i = 0; i < width; i++)
				principalComponents[j][i] += (left[i] - mean) * r;
		} // for j
	} // for s

	// multiply integer data
	for (int s = heightOfFloats; s < dim; s++)
	{
		Integer * left = static_cast<Integer *>(data[s]);
		Float mean = meanValsOfRows[s];

		for (int j = 0; j < numOfPrincipalComponents; j++)
		{
			Float r = featureVectors[s * numOfPrincipalComponents + j];// (s, j);

			for (int i = 0; i < width; i++)
				principalComponents[j][i] += (static_cast<Float>(left[i]) - mean) * r;
		} // for j
	} // for s
	return 0;
} // computePrincipalComponents

/* Principal Component Analysis of the data: naive approach, i.e. using correlation matrix.
     data - matrix that has float rows and integer rows, each row contains values of the same type:
            - data[0,...,heightOfFloats - 1] are pointers to float arrays;
            - data[heightOfFloats,...,heightOfFloats + heightOfInts - 1] are pointers to integer arrays;
     heightOfFloats - number of real-valued rows;
     heightOfInts - number of integer rows;
     width - width of each row;
	 numOfPrincipalComponents - number of principal components to be computed;
	 principalComponents - the principal components computed. */
int pca::pcaUsingCorrelationMatrix(void ** data,
	const int heightOfFloats,
	const int heightOfInts,
	const int width,
	const int numOfPrincipalComponents,
	Float ** principalComponents) noexcept
{
	int height = heightOfFloats + heightOfInts;

	// check number of principal components
	if (height < numOfPrincipalComponents)
		return 1;

	Vector<Float, Dynamic> meanValsOfRows(height);
	computeMeanOfEachRow(data, heightOfFloats, heightOfInts, width, meanValsOfRows.data());

	Matrix<Float, Dynamic, Dynamic> corMatrix(height, height);
	computeCorrelationMatrix(data, meanValsOfRows.data(), heightOfFloats, heightOfInts, width, corMatrix.data());
	
	// The following solver computes eigen vals & vectors: the order of eigen vals is increasing.
	SelfAdjointEigenSolver<Matrix<Float, Dynamic, Dynamic>> eigensolver(corMatrix);
 
	// Check result of eigen values & vectors computation.
	if (eigensolver.info() != Success)
		return 2; 

	// get feature vectors, taking into account increasing order of computed eigen values
	Matrix<Float, Dynamic, Dynamic, RowMajor> featureVectors
		= (eigensolver.eigenvectors().rowwise().reverse())(all, seq(0, numOfPrincipalComponents - 1));	

	computePrincipalComponents(data, meanValsOfRows.data(), featureVectors.data(),
		heightOfFloats, heightOfInts, width, numOfPrincipalComponents, principalComponents);

	return 0;
} // pcaNaive


// Compute mean value of array.
template<typename Type>
Float meanValue(Type * arr, const int length)
{
	Map<RowVector<Type, Dynamic>> vec(arr, length); // data assignment
	//return (Float)vec.mean(); //
	return static_cast<Float>(vec.sum()) / length;  // here, the method mean can be used, 
	                                                // but for the case of integers it provides integer mean!
} // meanValue

/* Compute mean value of each row of float and int Data.
     data - matrix that has float rows and integer rows, each row contains values of the same type:
            - data[0,...,heightOfFloats - 1] are pointers to float arrays;
            - data[heightOfFloats,...,heightOfFloats + heightOfInts - 1] are pointers to integer arrays;
     heightOfFloats - number of real-valued rows;
     heightOfInts - number of integer rows;
     width - width of each row;
     means - array of mean values that are computed. */
int pca::computeMeanOfEachRow(void ** data,
	const int heightOfFloats,
	const int heightOfInts,
	const int width,
	Float * means) noexcept
{
	// compute mean vals of real-valued rows using Eigen mean-func
	for (int i = 0; i < heightOfFloats; i++)
		means[i] = meanValue(static_cast<Float *>(data[i]), width);

	// compute mean vals of integer rows
	for (int i = heightOfFloats; i < heightOfFloats + heightOfInts; i++)
		means[i] = meanValue(static_cast<Integer *>(data[i]), width);

	return 0;
}

/* Compute matrix of correlations of the given data.
     data - matrix that has float rows and integer rows, each row contains values of the same type:
            - data[0,...,heightOfFloats - 1] are pointers to float arrays;
            - data[heightOfFloats,...,heightOfFloats + heightOfInts - 1] are pointers to integer arrays;
     means - array of mean values;
     heightOfFloats - number of real-valued rows;
     heightOfInts - number of integer rows;
     width - width of each row;
     correlations - matrix of the correlations computed. */
int pca::computeCorrelationMatrix(void ** data,
	const Float * means,
	const int heightOfFloats,
	const int heightOfInts,
	const int width,
	Float * correlations) noexcept
{
	/* Here, we use the formula:
	         C = 1/n * D * D.tr - mu.tr * mu,
       where C is d x d matrix of correlations computed,
	         D is an input d x n mtrix,
			 mu is a vector that contains mean values of the matrix D rows,
			 D.tr and mu.tr are transposed matrices.*/

	int dim = heightOfFloats + heightOfInts;
		
	for (int i = 0; i < heightOfFloats; i++)
	{
		Float * left = static_cast<Float *> (data[i]);

		// compute float-by-float correlations
		for (int j = i; j < heightOfFloats; j++)
		{
			Float * right = static_cast<Float *> (data[j]);
			Double res = 0.0;

			for (int k = 0; k < width; k++)
				res += left[k] * right[k];

			correlations[i * dim + j] = correlations[j * dim + i] 
				= static_cast<Float>( res / width - means[i] * means[j] );
		} // for j

		// compute float-by-integer correlations
		for (int j = 0; j < heightOfInts; j++)
		{
			Integer * right = static_cast<Integer *> (data[heightOfFloats + j]);
			Double res = 0.0;

			for (int k = 0; k < width; k++)
				res += left[k] * right[k];

			correlations[i * dim + heightOfFloats + j]
				= correlations[(heightOfFloats + j) * dim + i] 
				= static_cast<Float>( res / width - means[i] * means[heightOfFloats + j] );
		} // for j
	} // for i

	// compute integer-by-integer correlations
	for (int i = 0; i < heightOfInts; i++)
	{
		Integer * left = static_cast<Integer *> (data[heightOfFloats + i]);

		for (int j = i; j < heightOfInts; j++)
		{
			Integer res = 0;
			Integer * right = static_cast<Integer *> (data[heightOfFloats + j]);

			for (int k = 0; k < width; k++)
				res += left[k] * right[k];

			correlations[(i + heightOfFloats) * dim + heightOfFloats + j]
				= correlations[(j + heightOfFloats) * dim + heightOfFloats + i]
				= static_cast<Float>(res) / width - means[i + heightOfFloats] * means[j + heightOfFloats];
		} // for j
	} // for i

	return 0;
} // computeCorrelationMatrix

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
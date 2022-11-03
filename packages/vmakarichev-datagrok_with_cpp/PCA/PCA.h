// PCA.h
// Principal Component Analysis (PCA) using the lib Eigen: headers of functions

// REMARK 1. Each row of the input data contains Datagrok column. 
//           For this reason, the following convention is used:
//             -  height is a number of Datagrok columns to be processed,
//             -  width is a number of each Datagrok column.

// REMARK 2. Here, we operate matrices that have float rows and integer rows, 
//           each row contains values of the same type.
//           Each matrix consists of two blocks: float rows and integer rows.
//           In this case, an input is void **.

// RMEARK 3. Also, the same methods are implemented for the case when data is
//           given by float *.

#ifndef PCA_H
#define PCA_H

namespace pca {

	typedef float Float;
	typedef int Integer;
	typedef double Double;

	/* Principal Component Analysis of the data: using correlation matrix.
	     data - input matrix;
	     height, width - sizes of the input;
	     numOfPrincipalComponents - number of principal components to be computed;
	     principalComponents - the principal components computed;
		 approxData - approximation of the input data using principal components obtained.  */
	int pcaUsingCorrelationMatrix(Float * data,
		const int height,
		const int width,
		const int numOfPrincipalComponents,
		Float * principalComponents,
		Float * approxData = NULL) noexcept;

	/* Principal Component Analysis of the data: using correlation matrix.
	     data - matrix that has float rows and integer rows, each row contains values of the same type:
	            - data[0,...,heightOfFloats - 1] are pointers to float arrays;
	            - data[heightOfFloats,...,heightOfFloats + heightOfInts - 1] are pointers to integer arrays;
	     heightOfFloats - number of real-valued rows;
	     heightOfInts - number of integer rows;
	     width - width of each row;
	     numOfPrincipalComponents - number of principal components to be computed;
		 principalComponents - the principal components computed.  */
	int pcaUsingCorrelationMatrix(void ** data,		
		const int heightOfFloats,
		const int heightOfInts,
		const int width,
		const int numOfPrincipalComponents,
		Float ** principalComponents) noexcept;

	/* Compute mean value of each row of float and int Data.
	     data - matrix that has float rows and integer rows, each row contains values of the same type: 
		        - data[0,...,heightOfFloats - 1] are pointers to float arrays;
				- data[heightOfFloats,...,heightOfFloats + heightOfInts - 1] are pointers to integer arrays;
	     heightOfFloats - number of real-valued rows;
	     heightOfInts - number of integer rows;
	     width - width of each row;
	     means - array of mean values that are computed. */
	int computeMeanOfEachRow(void ** data,
		const int heightOfFloats,
		const int heightOfInts,
		const int width,
		Float * means) noexcept;
	
	/* Compute matrix of correlations of the given data.
	     data - matrix that has float rows and integer rows, each row contains values of the same type:
	            - data[0,...,heightOfFloats - 1] are pointers to float arrays;
	            - data[heightOfFloats,...,heightOfFloats + heightOfInts - 1] are pointers to integer arrays;
		 means - array of mean values;
	     heightOfFloats - number of real-valued rows;
	     heightOfInts - number of integer rows;
	     width - width of each row;
	     correlations - matrix of the correlations computed. */
	int computeCorrelationMatrix(void ** data,
		const Float * means,
		const int heightOfFloats,
		const int heightOfInts,
		const int width,
		Float * correlations) noexcept;

	// Maximum absolute deviation between arrays
	Float mad(Float * arr1, Float * arr2, const int length);
};

#endif // PCA_H


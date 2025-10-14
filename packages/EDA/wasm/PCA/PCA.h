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

	const Float TOL = 0.000001f;
	const Integer MAX_ITER = 100;

	enum ResultCode {NO_ERROR = 0, UNCORRECT_ARGUMENTS_ERROR, COMPUTATION_ERROR, METHOD_ERROR};

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
		const int centerNum,
		const int scaleNum,
		Float * principalComponents,
		Float * approxData = 0) noexcept;	

	// Maximum absolute deviation between arrays
	Float mad(Float * arr1, Float * arr2, const int length) noexcept;

    /* The NIPALS algorithm for PCA.
         data - input matrix;
         height, width - sizes of the input;
	     numOfPrincipalComponents - number of principal components to be computed;
	     principalComponents - the principal components computed.

	   Reference
	     H. Risvik. Principal component analysis (PCA) & NIPALS algorithm.
	     https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=a1a75889929b604cfefcb7cddfd58accce8d01d8
    */
	int nipals(Float* data,
	    const int height,
	    const int width,
	    const int numOfPrincipalComponents,
	    Float* principalComponents) noexcept;
};

#endif // PCA_H


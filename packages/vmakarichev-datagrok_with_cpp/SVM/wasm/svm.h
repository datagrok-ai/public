// svm.h

/* Implementation of the method SVM (Support Vector Machine) for DATAGROK.

   The following references are used:
   [1] Suykens, J., Vandewalle, J. "Least Squares Support Vector Machine Classifiers",
	   Neural Processing Letters 9, 293–300 (1999). https://doi.org/10.1023/A:1018628609742
*/

#ifndef SVM_H
#define SVM_H

#include <cmath>
using std::sqrt;
using std::fabs;

// TODO: this include should be removed
#include<iostream>
using namespace std;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

namespace svm {

	// computation result code
	enum ResultCode {
		NO_ERRORS = 0,
		UNKNOWN_PROBLEM,
		INCORRECT_HYPERPARAMETER,
		INCORRECT_PARAMETER_OF_KERNEL,
		UNKNOWN_KERNEL_TYPE,
		INCORRECT_KERNEL_PARAMS_COUNT,
		INCORRECT_PROBABILITY,
		INCORRECT_PERCENTAGE,
		INCORRECT_SIZE
	};

	// types of model kernels
	enum KernelType { LINEAR = 0, POLYNOMIAL, RBF, SIGMOID };

	const int MAX_NUM_OF_KERNEL_PARAM = 2;

	/* Check correctness of kernel parameters.
	      kernel - kernel type,
		  kernelParameters - parameters of kernel. */
	template<typename Float>
	bool areKernelParametersCorrect(int kernel, Float kernelParameters[MAX_NUM_OF_KERNEL_PARAM])
	{
		Float zero = static_cast<Float>(0);
		Float c = kernelParameters[0];
		Float d = kernelParameters[1];
		Float sigma = kernelParameters[0];

		switch (kernel)
		{
		case LINEAR:
			return true;
		case POLYNOMIAL:			
			return ((c > zero) && (d > zero));
		case RBF:			
			return (c > zero);
		case SIGMOID:
			return true;
		default:
			return false;
		}
	} // areKernelParametersCorrect

	/*  Value of SVM-kernel function.
		   kernel - type of kernel
		   kernelParams - parameters of kernel
		   v1, v2 - kernel arguments. */
	template<typename Float, typename VecType1, typename VecType2>
	Float K(int kernel, Float kernelParams[MAX_NUM_OF_KERNEL_PARAM], VecType1& v1, VecType2& v2)
	{
		switch (kernel)
		{
		case LINEAR:
			return v1.dot(v2);
		case RBF:
			return exp(-(v1 - v2).squaredNorm() / (kernelParams[0] * kernelParams[0]));
		}
		return 0;
	} // K

}; // svm

#endif // SVM_H



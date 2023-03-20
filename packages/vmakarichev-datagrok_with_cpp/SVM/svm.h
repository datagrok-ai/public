// svm.h

/* Implementation of the method SVM (Support Vector Machine) for DATAGROK.
   
   The following references are used:
   [1] Suykens, J., Vandewalle, J. "Least Squares Support Vector Machine Classifiers", 
	   Neural Processing Letters 9, 293–300 (1999). https://doi.org/10.1023/A:1018628609742
*/

#ifndef SVM_H
#define SVM_H

// TODO: this include should be removed
#include<iostream>
using namespace std;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

namespace svm {

	// computation result code
	enum ResultCode { NO_ERRORS = 0, 
		UNKNOWN_PROBLEM,
		INCORRECT_HYPERPARAMETER,
		INCORRECT_PARAMETER_OF_KERNEL,
		UNKNOWN_KERNEL_TYPE
	};

	// types of model kernels
	enum KernelType {LINEAR = 0, POLYNOMIAL, RBF, SIGMOID};

	const int MAX_NUM_OF_KERNEL_PARAM = 2;


	/* Create dataset from columns data.
	      columsDataPtr - pointer to columns data
		  rowCount - number of rows
		  colCount - number of columns
		  datasetPtr - pointer to dataset
		  
	   REMARK. In DATAGROK, column-oriented data storage is used,
	           but row-oriented approach is preffered in SVM, and
			   this function provides it.	*/
	template<typename Float>
	int createDataset(Float * columsDataPtr, int rowCount, int colCount, Float* datasetPtr) noexcept
	{	
		// pointers-to-matrices assignment
		Map < Matrix<Float, Dynamic, Dynamic, ColMajor>> A(columsDataPtr, rowCount, colCount);
		Map < Matrix<Float, Dynamic, Dynamic, RowMajor>> B(datasetPtr, rowCount, colCount);

		// get row-oriented data representation
		B = A;

		return NO_ERRORS;
	}


	/* Compute matrix of the linear system for the LS-SVM method 
	   with LINEAR kernel.
	      gammaInv - value inverse to hyperparameter gamma	      
		  xTrain - feature vectors for training model
		  yTrain - labels of feature vectors 
		  samplesCount - number of training samples
		  featuresCount - number of features, i.e. feature space dimension
		  A - matrix of linear system for the LS-SVM method  */
	template<typename Float>
	int computeMatrixOfLSSVMsystemWithLinearKernel(Float gammaInv,
		Float* xTrain, Float* yTrain, int samplesCount, int featuresCount,
		Matrix<Float, Dynamic, Dynamic, RowMajor>& A)
	{		
		// assign train data pointer with the matrix X
		Map<Matrix<Float, Dynamic, Dynamic, RowMajor>> X(xTrain, samplesCount, featuresCount);
		//cout << "\nX:\n" << X << endl;

		// compute left upper block
		for (int i = 0; i < samplesCount; i++)
		{
			for (int j = 0; j < i; j++)
				A(i, j) = A(j, i) = yTrain[i] * yTrain[j] * X.row(i).dot(X.row(j));

			A(i, i) = X.row(i).squaredNorm() + gammaInv;
		}
		// compute left lower and rigth upper block
		for (int j = 0; j < samplesCount; j++)
			A(samplesCount, j) = A(j, samplesCount) = yTrain[j];

		// right lower element
		A(samplesCount, samplesCount) = 0;

		return NO_ERRORS;
	} // computeMatrixOfLSSVMsystemWithLinearKernel


	/* Compute matrix of the linear system for the LS-SVM method.
	   Linear, polynomial, RBF and sigmoid kernels are considered.
	      gammaInv - value inverse to hyperparameter gamma
		  kernelType - type of the kernel applied
		  kernelParams - parameters of kernel
		  xTrain - feature vectors for training model
		  yTrain - labels of feature vectors 
		  samplesCount - number of training samples
		  featuresCount - number of features, i.e. feature space dimension
		  A - matrix of linear system for the LS-SVM method  */
	template<typename Float>
	int computeMatrixOfLSSVMsystem(Float gammaInv, int kernelType, Float* kernelParams,
		Float* xTrain, Float* yTrain, int samplesCount, int featuresCount,
		Matrix<Float, Dynamic, Dynamic, RowMajor>& A)
	{
		switch (kernelType)
		{
		case LINEAR:
			return computeMatrixOfLSSVMsystemWithLinearKernel(gammaInv, xTrain, yTrain, samplesCount, featuresCount, A);

		// TODO: consider other cases!
		default:
			return UNKNOWN_KERNEL_TYPE;
		} 
	} // computeMatrixOfLSSVMsystem


	/* Train Least Square Support Vector Machine (LS-SVM).
	      gamma - hyperparameter
	      kernelType - type of the kernel applied
		  kernelParams - parameters of kernel
		  xTrain - feature vectors for training model
		  yTrain - labels of feature vectors 
		  samplesCount - number of training samples
		  featuresCount - number of features, i.e. feature space dimension
		  modelParams - parameters of model that is trained  */
	template<typename Float>
	int trainLSSVM( Float gamma, int kernelType, Float kernelParams[MAX_NUM_OF_KERNEL_PARAM],
		Float* xTrain, Float* yTrain, int samplesCount, int featuresCount,
		Float* modelParams) noexcept
	{
		/* In order to find paramters of LS - SVM model, a special system of linear algebraic equations
		   should be solved (see [1] for more details). 
		   So, the following pricipal steps are:
		      1) compute the matrix of this system (A);
			  2) compute the rigth hand side of this system (b);
			  3) solve the system Ax = b. */

		// check gamma value
		if (gamma <= static_cast<Float>(0.0))
			return INCORRECT_HYPERPARAMETER;

		Float gammaInv = static_cast<Float>(1.0) / gamma;

		// matrix of the system to be further solved
		Matrix<Float, Dynamic, Dynamic, RowMajor> A(samplesCount + 1, samplesCount + 1);

		// compute the matrix A
		int resCode = computeMatrixOfLSSVMsystem(gammaInv, kernelType, kernelParams,
			xTrain, yTrain, samplesCount, featuresCount, A);

		// check results of the matrix A computation
		if (resCode != NO_ERRORS)
			return resCode;

		// create rigth hand side of 
		Vector<Float, Dynamic> b = Vector<Float, Dynamic>::Ones(samplesCount + 1);
		b(samplesCount) = 0;

		// assign modelParams with a vector
		Map<Vector<Float, Dynamic>> x(modelParams, samplesCount + 1);

		// solve the system required: 
		x = A.partialPivLu().solve(b);

		// The following is for testing!
		/*cout << "\nA:\n" << A << endl;

		cout << "\nrigth part:\n"
			<< b << endl;

		Vector<Float, Dynamic> params = A.partialPivLu().solve(b);
		
		cout << "\nsolution:\n" << params << endl;

		Vector<Float, Dynamic> w = Vector<Float, Dynamic>::Zero(featuresCount);

		Map<Matrix<Float, Dynamic, Dynamic, RowMajor>> X(xTrain, samplesCount, featuresCount);

		for (int i = 0; i < samplesCount; i++)
			w += params(i) * yTrain[i] * X.row(i);

		cout << "\nw:\n" << w << endl;*/

		return NO_ERRORS;
	} // trainModel

}; // svm

#endif // SVM_H


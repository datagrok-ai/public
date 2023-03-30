// svm.h

/* Implementation of the method SVM (Support Vector Machine) for DATAGROK.

   The following references are used:
   [1] Suykens, J., Vandewalle, J. "Least Squares Support Vector Machine Classifiers",
	   Neural Processing Letters 9, 293�300 (1999). https://doi.org/10.1023/A:1018628609742
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

	// constants for kernel params
	const int MAX_NUM_OF_KERNEL_PARAM = 2;
	const int RBF_SIGMA_INDEX = 0;
	const int POLYNOMIAL_C_INDEX = 0;
	const int POLYNOMIAL_D_INDEX = 1;
	const int SIGMOID_KAPPA_INDEX = 0;
	const int SIGMOID_THETA_INDEX = 1;

	// Check correctness of LS-SVM hyperparameter gamma	
	template<typename Float>
	bool isGammaCorrect(Float gamma) noexcept
	{
		return (gamma > static_cast<Float>(0));
	}

	/* Check correctness of kernel parameters.
	      kernel - kernel type,
		  kernelParameters - parameters of kernel. */
	template<typename Float>
	bool areKernelParametersCorrect(int kernel, Float kernelParameters[MAX_NUM_OF_KERNEL_PARAM]) noexcept
	{
		Float zero = static_cast<Float>(0);
		Float c = kernelParameters[POLYNOMIAL_C_INDEX];
		Float d = kernelParameters[POLYNOMIAL_D_INDEX];
		Float sigma = kernelParameters[RBF_SIGMA_INDEX];

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
	Float kernelFunc(int kernel, Float kernelParams[MAX_NUM_OF_KERNEL_PARAM], VecType1& v1, VecType2& v2) noexcept
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

	/* Compute matrix of the linear system for the LS-SVM method
	   with LINEAR kernel.
		  gammaInv - value inverse to hyperparameter gamma
		  xTrain - feature vectors for training model
		  yTrain - labels of feature vectors
		  samplesCount - number of training samples
		  featuresCount - number of features, i.e. feature space dimension
		  A - matrix of linear system for the LS-SVM method  */
	template<typename Float, typename MatrixType>
	int computeMatrixOfLSSVMsystemWithLinearKernel(Float gammaInv,
		Float* xTrain, Float* yTrain, int samplesCount, int featuresCount,
		MatrixType& A) noexcept
	{
		// assign train data pointer with the matrix X
		Map<Matrix<Float, Dynamic, Dynamic, RowMajor>> X(xTrain, samplesCount, featuresCount);

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

	/* Compute matrix of the linear system for the LS-SVM method
	   with RBF kernel.
		  gammaInv - value inverse to hyperparameter gamma
		  sigma - parameter of RBF kernel
		  xTrain - feature vectors for training model
		  yTrain - labels of feature vectors
		  samplesCount - number of training samples
		  featuresCount - number of features, i.e. feature space dimension
		  A - matrix of linear system for the LS-SVM method  */
	template<typename Float, typename MatrixType>
	int computeMatrixOfLSSVMsystemWithRBFkernel(Float gammaInv, Float sigma,
		Float* xTrain, Float* yTrain, int samplesCount, int featuresCount,
		MatrixType& A) noexcept
	{
		//cout << "\nsigma: " << sigma << endl;

		Float sigmaSquared = sigma * sigma;
		//cout << "\nsigma^2: " << sigmaSquared << endl;

		// assign train data pointer with the matrix X
		Map<Matrix<Float, Dynamic, Dynamic, RowMajor>> X(xTrain, samplesCount, featuresCount);
		//cout << "\nX:\n" << X << endl;

		// compute left upper block
		for (int i = 0; i < samplesCount; i++)
		{
			for (int j = 0; j < i; j++) {
				// this is for testing
				/*cout << "(" << i << ", " << j << "):\n"
					<< "   x_i: " << X.row(i) << endl
					<< "   x_j: " << X.row(j) << endl
					<< "   x_i - x_j: " << X.row(i) - X.row(j) << endl
					<< "   |x_i - x_j|^2: " << (X.row(i) - X.row(j)).squaredNorm() << endl
					<< "   -|x_i - x_j|^2 / sigma^2: " << -(X.row(i) - X.row(j)).squaredNorm() / sigmaSquared << endl
					<< "   exp(-|x_i - x_j|^2 / sigma^2): " << exp(-(X.row(i) - X.row(j)).squaredNorm() / sigmaSquared) << endl
					<< "   labels: " << yTrain[i] << ", " << yTrain[j] << endl
					<< "\n\n\n";*/

				A(i, j) = A(j, i) = yTrain[i] * yTrain[j] * exp(-(X.row(i) - X.row(j)).squaredNorm() / sigmaSquared);
			}				

			A(i, i) = gammaInv + 1; // here, 1 = exp(-|x_i - x_i|^2 / sigma^2) = exp( 0 )
		}
		// compute left lower and rigth upper block
		for (int j = 0; j < samplesCount; j++)
			A(samplesCount, j) = A(j, samplesCount) = yTrain[j];

		// right lower element
		A(samplesCount, samplesCount) = 0;

		// this is for testing
		//cout << "\nA:\n" << A << endl;

		return NO_ERRORS;
	} // computeMatrixOfLSSVMsystemWithRBFkernel

	/* Compute matrix of the linear system for the LS-SVM method
	   with POLYNOMIAL kernel.
		  gammaInv - value inverse to hyperparameter gamma
		  cParam - parameter of polynomial kernel (c)
		  dParam - parameter of polynomial kernel (d)
		  xTrain - feature vectors for training model
		  yTrain - labels of feature vectors
		  samplesCount - number of training samples
		  featuresCount - number of features, i.e. feature space dimension
		  A - matrix of linear system for the LS-SVM method  */
	template<typename Float, typename MatrixType>
	int computeMatrixOfLSSVMsystemWithPolynomialKernel(Float gammaInv, Float cParam, Float dParam,
		Float* xTrain, Float* yTrain, int samplesCount, int featuresCount,
		MatrixType& A) noexcept
	{
		/*cout << "\nc: " << cParam << endl;
		cout << "\nd: " << dParam << endl;*/		

		// assign train data pointer with the matrix X
		Map<Matrix<Float, Dynamic, Dynamic, RowMajor>> X(xTrain, samplesCount, featuresCount);
		//cout << "\nX:\n" << X << endl;

		// compute left upper block
		for (int i = 0; i < samplesCount; i++)
		{
			for (int j = 0; j < i; j++) {
				// this is for testing
				/*cout << "(" << i << ", " << j << "):\n"
					<< "   x_i: " << X.row(i) << endl
					<< "   x_j: " << X.row(j) << endl
					<< "   <x_i,  x_j>: " << X.row(i).dot(X.row(j)) << endl
					<< "   1 + <x_i,  x_j> / c: " << X.row(i).dot(X.row(j)) / cParam + 1 << endl
					<< "   (1 + <x_i,  x_j> / c)^d: " << pow(X.row(i).dot(X.row(j)) / cParam + 1, dParam) << endl					
					<< "   labels: " << yTrain[i] << ", " << yTrain[j] << endl
					<< "\n\n\n";*/

				//  here, 1 is used by the formula from [1]
				A(i, j) = A(j, i) = yTrain[i] * yTrain[j] * pow(X.row(i).dot(X.row(j)) / cParam + 1, dParam);
			}

			//  here, 1 is used by the formula from [1]
			A(i, i) = gammaInv + pow(X.row(i).squaredNorm() / cParam + 1, dParam);
		}
		// compute left lower and rigth upper block
		for (int j = 0; j < samplesCount; j++)
			A(samplesCount, j) = A(j, samplesCount) = yTrain[j];

		// right lower element
		A(samplesCount, samplesCount) = 0;

		// this is for testing
		//cout << "\nA:\n" << A << endl;

		return NO_ERRORS;
	} // computeMatrixOfLSSVMsystemWithPolynomialKernel

	/* Compute matrix of the linear system for the LS-SVM method
	   with SIGMOID kernel.
		  gammaInv - value inverse to hyperparameter gamma
		  cParam - parameter of polynomial kernel (c)
		  dParam - parameter of polynomial kernel (d)
		  xTrain - feature vectors for training model
		  yTrain - labels of feature vectors
		  samplesCount - number of training samples
		  featuresCount - number of features, i.e. feature space dimension
		  A - matrix of linear system for the LS-SVM method  */
	template<typename Float, typename MatrixType>
	int computeMatrixOfLSSVMsystemWithSigmoidKernel(Float gammaInv, Float kappa, Float theta,
		Float* xTrain, Float* yTrain, int samplesCount, int featuresCount,
		MatrixType& A) noexcept
	{
		/*cout << "\nc: " << cParam << endl;
		cout << "\nd: " << dParam << endl;*/

		// assign train data pointer with the matrix X
		Map<Matrix<Float, Dynamic, Dynamic, RowMajor>> X(xTrain, samplesCount, featuresCount);
		//cout << "\nX:\n" << X << endl;

		// compute left upper block
		for (int i = 0; i < samplesCount; i++)
		{
			for (int j = 0; j < i; j++) {
				// this is for testing
				/*cout << "(" << i << ", " << j << "):\n"
					<< "   x_i: " << X.row(i) << endl
					<< "   x_j: " << X.row(j) << endl
					<< "   <x_i,  x_j>: " << X.row(i).dot(X.row(j)) << endl
					<< "   1 + <x_i,  x_j> / c: " << X.row(i).dot(X.row(j)) / cParam + 1 << endl
					<< "   (1 + <x_i,  x_j> / c)^d: " << pow(X.row(i).dot(X.row(j)) / cParam + 1, dParam) << endl
					<< "   labels: " << yTrain[i] << ", " << yTrain[j] << endl
					<< "\n\n\n";*/

					//  here, 1 is used by the formula from [1]
				A(i, j) = A(j, i) = yTrain[i] * yTrain[j] * tanh(kappa * X.row(i).dot(X.row(j)) + theta);
			}

			//  here, 1 is used by the formula from [1]
			A(i, i) = gammaInv + tanh(kappa * X.row(i).squaredNorm() + theta);
		}
		// compute left lower and rigth upper block
		for (int j = 0; j < samplesCount; j++)
			A(samplesCount, j) = A(j, samplesCount) = yTrain[j];

		// right lower element
		A(samplesCount, samplesCount) = 0;

		// this is for testing
		//cout << "\nA:\n" << A << endl;

		return NO_ERRORS;
	} // computeMatrixOfLSSVMsystemWithSigmoidKernel

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
	template<typename Float, typename MatrixType>
	int computeMatrixOfLSSVMsystem(Float gammaInv, int kernelType, Float* kernelParams,
		Float* xTrain, Float* yTrain, int samplesCount, int featuresCount,
		MatrixType& A) noexcept
	{
		switch (kernelType)
		{
		case LINEAR: // linear kernel case
			return computeMatrixOfLSSVMsystemWithLinearKernel(gammaInv, xTrain, yTrain, samplesCount, featuresCount, A);

		case RBF: // RBF kernel case
			return computeMatrixOfLSSVMsystemWithRBFkernel(gammaInv, kernelParams[RBF_SIGMA_INDEX],
				xTrain, yTrain, samplesCount, featuresCount, A);

		case POLYNOMIAL: // polynomial kernel case
			return computeMatrixOfLSSVMsystemWithPolynomialKernel(gammaInv,
				kernelParams[POLYNOMIAL_C_INDEX], kernelParams[POLYNOMIAL_D_INDEX],
				xTrain, yTrain, samplesCount, featuresCount, A);

		case SIGMOID: // sigmoid kernel case
			return computeMatrixOfLSSVMsystemWithSigmoidKernel(gammaInv,
				kernelParams[SIGMOID_KAPPA_INDEX], kernelParams[SIGMOID_THETA_INDEX],
				xTrain, yTrain, samplesCount, featuresCount, A);

		default:
			return UNKNOWN_KERNEL_TYPE;
		}
	} // computeMatrixOfLSSVMsystem

	/* Train Least Square Support Vector Machine (LS-SVM).
		  gamma - hyperparameter
		  kernelType - type of the kernel applied
		  kernelParams - parameters of kernel
		  xTrain - feature vectors for training model (normalized)
		  yTrain - labels of feature vectors
		  samplesCount - number of training samples
		  featuresCount - number of features, i.e. feature space dimension
		  modelParams - parameters of model that is trained
		  weights - weights of the cos function, they computed in the case of linear kernel

	   WARNING. Training data should be normalized, i.e. mean value of each feature is 0
				and standard deviantion is 1.
	*/
	template<typename Float>
	int trainLSSVM(Float gamma, int kernel, Float kernelParams[MAX_NUM_OF_KERNEL_PARAM],
		Float* xTrain, Float* yTrain, int samplesCount, int featuresCount,
		Float* modelParams, Float* weights) noexcept
	{
		/* In order to find paramters of LS - SVM model, a special system of linear algebraic equations
		   should be solved (see [1] for more details).
		   So, the following pricipal steps are:
			  1) compute the matrix of this system (A);
			  2) compute the rigth hand side of this system (b);
			  3) solve the system Ax = b.

		   Also, model weights are computed in the case of linear kernel. */

		// check gamma value
		if (gamma <= static_cast<Float>(0.0))
			return INCORRECT_HYPERPARAMETER;

		//check kernel parameters
		if (!areKernelParametersCorrect(kernel, kernelParams))
			return INCORRECT_PARAMETER_OF_KERNEL;

		Float gammaInv = static_cast<Float>(1.0) / gamma;

		// matrix of the system to be further solved
		Matrix<Float, Dynamic, Dynamic, RowMajor> A(samplesCount + 1, samplesCount + 1);

		// compute the matrix A
		int resCode = computeMatrixOfLSSVMsystem(gammaInv, kernel, kernelParams,
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
		x = A.fullPivLu().solve(b);

		// The following is for testing!
		/*cout << "\nsolution:\n" << x << endl;
		cout << "\ndeviation:\n" << (A * x - b).norm() << endl;*/

		// finish computations in the case of non-linear kernel
		if (kernel != LINEAR)
			return NO_ERRORS;

		// compute bias 
		weights[featuresCount] = modelParams[samplesCount];

		// assign weights with w-vector
		Map<Vector<Float, Dynamic>> w(weights, featuresCount);

		// initialization 
		w = Vector<Float, Dynamic>::Zero(featuresCount);

		// assign normalized train data with 
		Map<Matrix<Float, Dynamic, Dynamic, RowMajor>> X(xTrain, samplesCount, featuresCount);

		// compute weigths
		for (int i = 0; i < samplesCount; i++)
			w += x(i) * yTrain[i] * X.row(i);

		// The following is for testing!
		/*cout << "\nsolution:\n" << x << endl;
		cout << "\nw:\n" << w << endl;*/

		return NO_ERRORS;
	} // trainModel

	/*  Predict labels of the target data using linear kernel model and precomputed weigths.		  
		  precomputedWeights - precomputed weights of the model
		  targetDataMatrix - matrix of the target data
		  prediction - target labels
		  targetSamplesCount - number of target samples	*/
	template<typename Float, typename MatrixType>
	int predictByLSSVMwithLinearKernel(Float* precomputedWeights,
		MatrixType & targetDataMatrix, Float* prediction) noexcept
	{	
		// get target data sizes		
		auto targetSamplesCount = targetDataMatrix.rows();
		auto featuresCount = targetDataMatrix.cols();

		// bias of the model
		Float bias = precomputedWeights[featuresCount];

		// assign weights-vector with the corresponding data pointer
		Map<Vector<Float, Dynamic>> w(precomputedWeights, featuresCount);		

		// predict labels of the target data
		for (int i = 0; i < targetSamplesCount; i++)
		{
			// put target data to the model
			Float cost = bias + w.dot(targetDataMatrix.row(i));

			// check sign and get label (see [1] for more details)
			if (cost > static_cast<Float>(0))
				prediction[i] = static_cast<Float>(1);
			else
				prediction[i] = static_cast<Float>(-1);
		}

		return NO_ERRORS;
	} // predictByLSSVMwithLinearKernel
		
	/**/
	template<typename Float, typename MatrixType>
	int predictByLSSVMwithRBFkernel(Float sigma, 
		Float* xTrain, Float* yTrain, int trainSamplesCount, int featuresCount, 
		Float* modelParams, MatrixType& targetDataMatrix, Float* prediction)
	{
		//cout << "\nsigma: " << sigma << endl;

		Float sigmaSquared = sigma * sigma;

		//cout << "\nsigma^2: " << sigmaSquared << endl;

		// assign normalized traine data matrix with normalized train data pointer
		Map<Matrix<Float, Dynamic, Dynamic, RowMajor>> X(xTrain, trainSamplesCount, featuresCount);

		/*cout << "\nNTD:\n" << targetDataMatrix << endl;
		cout << "\nX:\n" << X << endl;*/

		auto targetSamplesCount = targetDataMatrix.rows();

		//cout << "\ntarget samples count: " << targetSamplesCount << endl;

		// bias of the model
		Float bias = modelParams[trainSamplesCount];
		//cout << "\nbias: " << bias << endl;

		for (int j = 0; j < targetSamplesCount; j++)
		{
			// put target data to the model
			Float cost = bias;

			for(int i = 0; i < trainSamplesCount; i++)
				cost += modelParams[i] * yTrain[i] * exp(-(X.row(i) - targetDataMatrix.row(j)).squaredNorm() / sigmaSquared);

			// check sign and get label (see [1] for more details)
			if (cost > static_cast<Float>(0))
				prediction[j] = static_cast<Float>(1);
			else
				prediction[j] = static_cast<Float>(-1);
		} // for j

		return NO_ERRORS;
	} // predictByLSSVMwithRBFkernel

	/**/
	template<typename Float, typename MatrixType>
	int predictByLSSVMwithPolynomialKernel(Float cParam, Float dParam,
		Float* xTrain, Float* yTrain, int trainSamplesCount, int featuresCount,
		Float* modelParams, MatrixType& targetDataMatrix, Float* prediction)
	{
		/*cout << "\nc: " << cParam << endl;
		cout << "\nd: " << dParam << endl;*/

		// assign normalized traine data matrix with normalized train data pointer
		Map<Matrix<Float, Dynamic, Dynamic, RowMajor>> X(xTrain, trainSamplesCount, featuresCount);

		/*cout << "\nNTD:\n" << targetDataMatrix << endl;
		cout << "\nX:\n" << X << endl;*/

		auto targetSamplesCount = targetDataMatrix.rows();

		//cout << "\ntarget samples count: " << targetSamplesCount << endl;

		// bias of the model
		Float bias = modelParams[trainSamplesCount];
		//cout << "\nbias: " << bias << endl;

		for (int j = 0; j < targetSamplesCount; j++)
		{
			// put target data to the model
			Float cost = bias;

			for (int i = 0; i < trainSamplesCount; i++)				
				cost += modelParams[i] * yTrain[i] * pow(X.row(i).dot(targetDataMatrix.row(j)) / cParam + 1, dParam);				

			// check sign and get label (see [1] for more details)
			if (cost > static_cast<Float>(0))
				prediction[j] = static_cast<Float>(1);
			else
				prediction[j] = static_cast<Float>(-1);
		} // for j

		return NO_ERRORS;
	} // predictByLSSVMwithPolynomialKernel

	/**/
	template<typename Float, typename MatrixType>
	int predictByLSSVMwithSigmoidKernel(Float kappa, Float theta,
		Float* xTrain, Float* yTrain, int trainSamplesCount, int featuresCount,
		Float* modelParams, MatrixType& targetDataMatrix, Float* prediction)
	{
		// assign normalized traine data matrix with normalized train data pointer
		Map<Matrix<Float, Dynamic, Dynamic, RowMajor>> X(xTrain, trainSamplesCount, featuresCount);

		/*cout << "\nNTD:\n" << targetDataMatrix << endl;
		cout << "\nX:\n" << X << endl;*/

		auto targetSamplesCount = targetDataMatrix.rows();

		//cout << "\ntarget samples count: " << targetSamplesCount << endl;

		// bias of the model
		Float bias = modelParams[trainSamplesCount];
		//cout << "\nbias: " << bias << endl;

		for (int j = 0; j < targetSamplesCount; j++)
		{
			// put target data to the model
			Float cost = bias;

			for (int i = 0; i < trainSamplesCount; i++)
				cost += modelParams[i] * yTrain[i] * tanh(kappa * X.row(i).dot(targetDataMatrix.row(j)) + theta);				

			// check sign and get label (see [1] for more details)
			if (cost > static_cast<Float>(0))
				prediction[j] = static_cast<Float>(1);
			else
				prediction[j] = static_cast<Float>(-1);
		} // for j

		return NO_ERRORS;
	} // predictByLSSVMwithSigmoidKernel

	/* Predict labels using LS-SVM model.
		  kernelType - type of the kernel
		  kernelParams - parameters of kernel
		  xTrain - feature vectors for training model (normalized)
		  yTrain - labels of training vectors
		  trainSamplesCount - number of training vectors
		  featuresCount - number of features
		  means - mean values of training data features
		  stdDevs - standard deviations of training data features
		  modelParams - parameters of models
		  precomputedWeights - weights of the hyperplane
		  targetData - target data
		  labels - target labels
		  targetSamplesCount - number of target samples

	   REMARK. Precomputed weights reduce computations in the case of the linear kernel.
			   In other cases, model parameters (modelParams) are used.	*/
	template<typename Float>
	int predictByLSSVM(int kernelType, Float kernelParams[MAX_NUM_OF_KERNEL_PARAM],
		Float* xTrain, Float* yTrain, int trainSamplesCount, int featuresCount,
		Float* means, Float* stdDevs, Float* modelParams, Float* precomputedWeights,
		Float* targetData, Float* prediction, int targetSamplesCount) noexcept
	{
		//check kernel parameters
		if (!areKernelParametersCorrect(kernelType, kernelParams))
			return INCORRECT_PARAMETER_OF_KERNEL;

		// assign target data matrix (TD) with target data pointer
		Map<Matrix<Float, Dynamic, Dynamic, ColMajor>> TD(targetData, targetSamplesCount, featuresCount);

		// assign means with a vector
		Map<Vector<Float, Dynamic>> mu(means, featuresCount);

		// matrix for normalized target data (NTD)
		Matrix<Float, Dynamic, Dynamic, RowMajor> NTD(targetSamplesCount, featuresCount);

		// center target data
		NTD = TD.rowwise() - mu.transpose();

		// normilize target data
		for (int i = 0; i < featuresCount; i++)
		{
			Float sigma = stdDevs[i];

			if (sigma > static_cast<Float>(0))
				NTD.col(i) /= sigma;
		}		

		// compute prediction
		switch (kernelType)
		{
		case LINEAR: // the case of linear kernel
			return predictByLSSVMwithLinearKernel(precomputedWeights, NTD, prediction);

		case RBF: // the case of RBF kernel
			return predictByLSSVMwithRBFkernel(kernelParams[RBF_SIGMA_INDEX],
				xTrain, yTrain, trainSamplesCount, featuresCount, modelParams,
				NTD, prediction);

		case POLYNOMIAL: // the case of polynomial kernel
			return predictByLSSVMwithPolynomialKernel(
				kernelParams[POLYNOMIAL_C_INDEX], kernelParams[POLYNOMIAL_D_INDEX],
				xTrain, yTrain, trainSamplesCount, featuresCount, modelParams,
				NTD, prediction);

		case SIGMOID:
			return predictByLSSVMwithSigmoidKernel(
				kernelParams[SIGMOID_KAPPA_INDEX], kernelParams[SIGMOID_THETA_INDEX],
				xTrain, yTrain, trainSamplesCount, featuresCount, modelParams,
				NTD, prediction);

		default:
			return UNKNOWN_KERNEL_TYPE;
		}
	} // predictByLSSVM
}; // svm

#endif // SVM_H



// dataGenerators.h

// Tools for generating datasets for testing SVM.

#ifndef DATA_GENERATORS_H
#define DATA_GENERATORS_H

#include<cstdlib>
using namespace std;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "svm.h"

namespace svm
{
	// Constants for random generation
	const unsigned SEED = 10214313;
	const int RAND_SCALE = 1000;

	/* Change data labels by opposite values.
	   Each label value is replaced by the corresponding opposite one
	   with the specified probability
		  labels - data labels
		  samplesCount - number of labels
		  changeProbability - probability that each label is changed */
	template<typename Float>
	int changeLabels(Float* labels, int samplesCount, Float changeProbability) noexcept
	{
		using namespace svm;

		// check probability value
		if ((changeProbability < static_cast<Float>(0)) ||
			(changeProbability > static_cast<Float>(1)))
			return INCORRECT_PROBABILITY;

		// check size
		if (samplesCount < 1)
			return INCORRECT_SIZE;

		// randomize
		srand(SEED + samplesCount);

		// change values in a random manner
		for (int i = 0; i < samplesCount; i++)
			if (static_cast<Float>(rand() % RAND_SCALE) / RAND_SCALE < changeProbability)			
				labels[i] = -labels[i];
		
		return NO_ERRORS;
	} // changeLabels

	/* Generate dataset: separable case. Features are generated randomly using the uniform distribution.
       Each feature belongs to the corresponding segment [min, max].
	      kernel - type of kernel
	      kernelParams - parameters of kernel
	      featuresCount - number of features, i.e. dimension
	      samplesCount - number of the generated samples
	      minVal - min value
	      maxVal - max value
	      data - generated data
	      labels - generated labels
    WARNING. Memory for data and labels must be allocated outside this function.  */
	template<typename Float>
	int generateSeparable(int kernel, float kernelParams[MAX_NUM_OF_KERNEL_PARAM],
		int featuresCount, int samplesCount,
		Float minVal, Float maxVal,
		Float* data, Float* labels) noexcept
	{
		using namespace svm;

		// check parameters correctness
		if (!areKernelParametersCorrect(kernel, kernelParams))
			return INCORRECT_PARAMETER_OF_KERNEL;

		// check sizes
		if ((featuresCount < 1) || (samplesCount < 1))
			return INCORRECT_SIZE;

		// randomize
		srand(SEED + samplesCount + featuresCount);

		// assign data pointer with a matrix
		Map<Matrix<Float, Dynamic, Dynamic, ColMajor>> X(data, samplesCount, featuresCount);

		// generate random matrix: values from [-1, 1] are generated
		X = Matrix<Float, Dynamic, Dynamic, ColMajor>::Random(samplesCount, featuresCount);

		// generate core vector
		RowVector<float, Dynamic> v(featuresCount);		

		// linear transform coefficients
		Float c1 = (maxVal - minVal) / 2;
		Float c2 = (maxVal + minVal) / 2;

		// rescale data: each feature should belong to the correspondent [min, max] segment
		for (int i = 0; i < featuresCount; i++)
		{
			// linear [-1,1]-to-[min,max] transform
			X.col(i) = X.col(i) * c1 + c2 * Vector<Float, Dynamic>::Ones(samplesCount);

			Float randNum = static_cast<Float>(-0.5) + static_cast<Float>(rand() % RAND_SCALE) / RAND_SCALE;
			v(i) = randNum * c1 + c2;
		}

		// bias value
		Float bias = kernelFunc(kernel, kernelParams, v, v);

		// This is a heruistics
		if (kernel == RBF)
			bias /= 2;

		// auxilliry vector
		RowVector<float, Dynamic> w(featuresCount);

		// compute labels
		for (int i = 0; i < samplesCount; i++)
		{
			w = X.row(i);
			Float val = kernelFunc(kernel, kernelParams, w, v) - bias;

			labels[i] = (val > static_cast<Float>(0)) ? static_cast<Float>(1) : static_cast<Float>(-1);
		}

		return NO_ERRORS;
	} // generateSeparable

	/* Generate dataset: non-separable case.
	   Features are generated randomly using the uniform distribution.
	   Each feature belongs to the corresponding segment [min, max].
		  kernel - type of kernel
		  kernelParams - parameters of kernel
		  featuresCount - number of features, i.e. dimension
		  samplesCount - number of the generated samples
		  minVal - min value
		  maxVal - max value
		  data - generated data
		  labels - generated labels
		  violatorsPercentage - percentage of values that violate separability

	   WARNINGS. 1. Memory for data and labels must be allocated outside this function.
				 2. Since violators are generated randomly, actual number of vilators
					may differ from the given percentage. */
	template<typename Float>
	int generateNonSeparable(int kernel, float kernelParams[MAX_NUM_OF_KERNEL_PARAM],
		int featuresCount, int samplesCount,
		Float minVal, Float maxVal,
		Float* data, Float* labels,
		Float violatorsPercentage) noexcept
	{
		using namespace svm;

		// check percentage
		if ((violatorsPercentage < static_cast<Float>(0)) ||
			(violatorsPercentage > static_cast<Float>(100)))
			return INCORRECT_PERCENTAGE;

		// generate separable dataset
		int resCode = generateSeparable(kernel, kernelParams, featuresCount, samplesCount,
			minVal, maxVal, data, labels);
		if (resCode != NO_ERRORS)
			return resCode;

		// create violators
		return changeLabels(labels, samplesCount, violatorsPercentage / 100);
	} // generateNonSeparable
}; // svm

#endif // DATA_GENERATORS_H

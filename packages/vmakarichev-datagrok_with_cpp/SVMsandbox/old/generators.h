// generators.h

// Tools for generating test datasets for SVM

#include<iostream>
#include<cstdlib>
using namespace std;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

namespace gener
{
	// Result codes of generating procedures
	enum GenCode {NO_ERRORS = 0, 
		UNKNOWN_PROBLEM,
		INCORRECT_PROBABILITY,
		INCORRECT_PERCENTAGE
	};

	// Constants for random generation
	const unsigned SEED = 10214313;
	const int RAND_SCALE = 1000;

	/* Generate dataset: vectors are linearly separable by the hyperplane x * w + b = 0.
	   Features are generated randomly using the uniform distribution.
	   Each feature belongs to the corresponding segment [min, max].
	      weights - coordinates of the norm vector w
		  bias - value of b
		  featuresCount - number of features, i.e. dimension
		  samplesCount - number of the generated samples
		  minVals - array of min's
		  maxVals - array of max's
		  data - generated data
		  labels - generated labels

	   WARNING. Memory for data and labels must be allocated outside this function.  */
	template<typename Float>
	int generateLinearSeparable(Float * weights, Float bias, int featuresCount, int samplesCount, 
		Float * minVals, Float * maxVals,
		Float * data, Float * labels)
	{
		// TODO: add check of inputs correctness

		// assign data pointer with a matrix
		Map<Matrix<Float, Dynamic, Dynamic, ColMajor>> X(data, samplesCount, featuresCount);

		// generate random matrix: values from [-1, 1] are generated
		X = Matrix<Float, Dynamic, Dynamic, ColMajor>::Random(samplesCount, featuresCount);

		//cout << "\ngenerated X:\n" << X << endl;

		// rescale data: each feature should belong to the correspondent [min, max] segment
		for (int i = 0; i < featuresCount; i++)
		{
			Float c1 = (maxVals[i] - minVals[i]) / 2;
			Float c2 = (maxVals[i] + minVals[i]) / 2;

			// linear [-1,1]-to-[min,max] transform
			X.col(i) = X.col(i) * c1 + c2 * Vector<Float, Dynamic>::Ones(samplesCount);
		}

		//cout << "\nX after rescaling:\n" << X << endl;

		// assign weights pointer with a matrix 
		Map<Matrix<Float, Dynamic, Dynamic>> W(weights, featuresCount, 1);

		//cout << "\nW:\n" << W << endl;

		// put randomly generated points to the given hyperplane (x * w + b = 0) 
		Vector<Float, Dynamic> P = X * W + bias * Vector<Float, Dynamic>::Ones(samplesCount);

		//cout << "\nP:\n" << P << endl;

		// compute labels
		for (int i = 0; i < samplesCount; i++)
			labels[i] = (P(i) > static_cast<Float>(0)) ? static_cast<Float>(1) : static_cast<Float>(-1);

		return NO_ERRORS;
	} // generateLinearSeparable

	/* Change data labels by opposite values.
	   Each label value is replaced by the corresponding opposite one 
	   with the specified probability
	      labels - data labels
		  samplesCount - number of labels
		  changeProbability - probability that each label is changed  
		  changesPercentage - actual changes percentage */
	template<typename Float>
	int changeLabels(Float* labels, int samplesCount, 
		Float changeProbability, Float & changesPercentage)
	{
		// check probability value
		if ((changeProbability < static_cast<Float>(0)) ||
			(changeProbability > static_cast<Float>(1)))
			return INCORRECT_PROBABILITY;

		srand(SEED);

		int changesCount = 0;

		// change values in a random manner
		for (int i = 0; i < samplesCount; i++)
			if (static_cast<Float>(rand() % RAND_SCALE) / RAND_SCALE < changeProbability)
			{
				labels[i] = -labels[i];
				changesCount++;
			}

		// compute actual percentage of changes
		changesPercentage = static_cast<Float>(changesCount) / samplesCount * 100;

		return NO_ERRORS;
	} // changeLabels

	/* Generate dataset: vectors are not linearly separable by the hyperplane x * w + b = 0.
	   Features are generated randomly using the uniform distribution.
	   Each feature belongs to the corresponding segment [min, max].
		  weights - coordinates of the norm vector w
		  bias - value of b
		  featuresCount - number of features, i.e. dimension
		  samplesCount - number of the generated samples
		  minVals - array of min's
		  maxVals - array of max's
		  data - generated data
		  labels - generated labels
		  violatorsPercentage - percentage of values that violate separability
		  percentageOfActualViolators -actual percentage of violators

	   WARNINGS. 1. Memory for data and labels must be allocated outside this function.  
	             2. Since violators are generated randomly, actual number of vilators
				    may differ from the given percentage. */
	template<typename Float>
	int generateLinearNonSeparable(Float* weights, Float bias, int featuresCount, int samplesCount,
		Float* minVals, Float* maxVals,
		Float* data, Float* labels, 
		Float violatorsPercentage, Float & percentageOfActualViolators)
	{
		// check percentage
		if ((violatorsPercentage < static_cast<Float>(0)) ||
			(violatorsPercentage > static_cast<Float>(100)))
			return INCORRECT_PERCENTAGE;

		// generate linearly separable dataset
		int resCode = generateLinearSeparable(weights, bias, featuresCount, samplesCount,
			minVals, maxVals, data, labels);
		if (resCode != NO_ERRORS)
			return resCode;

		// create violators
		return changeLabels(labels, samplesCount, violatorsPercentage / 100, percentageOfActualViolators);
	} // generateLinearNonSeparable

	/**/
	template<typename Float>
	int generateModelParams(Float* weights, Float& bias, int featuresCount)
	{
		srand(SEED);

		for (int i = 0; i < featuresCount; i++)
			weights[i] = static_cast<Float>(rand() % (2 * RAND_SCALE) - RAND_SCALE) / RAND_SCALE;

		bias = static_cast<Float>(rand() % (2 * RAND_SCALE) - RAND_SCALE) / RAND_SCALE;

		return NO_ERRORS;
	} // generateModelParams

}; // gener

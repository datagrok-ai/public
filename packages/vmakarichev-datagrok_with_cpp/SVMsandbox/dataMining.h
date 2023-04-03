// dataMining.h

// Data mining tools

#ifndef DATA_MINING_H
#define DATA_MINING_H

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

// data mining tools
namespace dmt {

	enum ResultCode { 
		NO_ERRORS = 0, 
		UNKNOWN_PROBLEM,
		INCORRECT_SIZE
   	};

	// confusion matrix constants
	const int CONFUSION_MATR_SIZE = 4;
	const int TRUE_POSITIVE_INDEX = 0;
	const int FALSE_NEGATIVE_INDEX = 1;
	const int FALSE_POSITIVE_INDEX = 2;
	const int TRUE_NEGATIVE_INDEX = 3;

	/* Create normalized dataset from columns data.
	   Each column of the ouput is centered and normalized.
		  columsData - pointer to columns data
		  rowCount - number of rows
		  colCount - number of columns
		  normalizedDataRows - pointer to normalized data rows
		  means - mean values of source columns
		  stdDevs - standard deviations of source columns

	 REMARKS. 1. In DATAGROK, column-oriented data storage is used,
	        	 but row-oriented approach is preffered in SVM, and
				 this function provides it.
			  2. Row-oriented data storage is a result. */
	template<typename Float>
	int getNormalizedDataset(Float* columsData, int rowCount, int colCount,
		Float* normalizedDataRows, Float* means, Float* stdDevs) noexcept
	{
		// check sizes
		if ((rowCount < 1) || (colCount < 1))
			return INCORRECT_SIZE;

		// pointers-to-matrices assignment
		Map < Matrix<Float, Dynamic, Dynamic, ColMajor>> A(columsData, rowCount, colCount);
		Map < Matrix<Float, Dynamic, Dynamic, RowMajor>> B(normalizedDataRows, rowCount, colCount);
		Map < Vector<Float, Dynamic> > mu(means, colCount);
		Map < Vector<Float, Dynamic> > sigma(stdDevs, colCount);

		// compute mean values & standard deviations
		for (int i = 0; i < colCount; i++)
		{
			mu(i) = A.col(i).mean();
			sigma(i) = sqrt(A.col(i).squaredNorm() / rowCount - mu(i) * mu(i));
		}

		// get A centered
		B = A.rowwise() - mu.transpose();

		// norm columns of B
		for (int i = 0; i < colCount; i++)
		{
			Float current = sigma(i);

			if (current > static_cast<Float>(0))
				B.col(i) /= current;
		}		

		return NO_ERRORS;
	} // createNormalizedDataset

	/* Compare labels and their prediciotns: BINARY CLASSIFICATION CASE.
	      labels - training labels
		  predictions - predicted labels
		  correctness - array of mistakes (1 - correct prediction, 0 - incorrect prediction)
		  samplesCount - number of training samples
		  confusionMatrix - confusion matrix  */
	template<typename Float>
	int compareLabelsAndTheirPredictions(Float* labels, Float* predictions,
		Float* correctness, int samplesCount,
		int confusionMatrix[CONFUSION_MATR_SIZE])
	{
		Float zero = static_cast<Float>(0);

		// initialization
		for (int i = 0; i < CONFUSION_MATR_SIZE; i++)
			confusionMatrix[i] = 0;

		// labels vs. prediction comparison 
		for (int i = 0; i < samplesCount; i++)
		{
			correctness[i] = labels[i] * predictions[i];

			if (labels[i] > zero)
				if (predictions[i] > zero)
					confusionMatrix[TRUE_POSITIVE_INDEX]++;
				else
					confusionMatrix[FALSE_NEGATIVE_INDEX]++;
			else
				if (predictions[i] > zero)
					confusionMatrix[FALSE_POSITIVE_INDEX]++;
				else
					confusionMatrix[TRUE_NEGATIVE_INDEX]++;
		}

		return NO_ERRORS;
	} // compareLabelsAndTheirPredictions

} // dmt

#endif // DATA_MINING_H


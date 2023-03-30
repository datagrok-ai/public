// dataMining.h

// Data mining tools

#ifndef DATA_MINING_H
#define DATA_MINING_H

#include<iostream>
using namespace std;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

// data mining tools
namespace dmt {

	enum ResultCode { 
		NO_ERRORS = 0, 
		UNKNOWN_PROBLEM,
		INCORRECT_SIZE
   	};

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

		// For testing
		/*cout << "\nA:\n" << A
			<< "\n\nB:\n" << B
			<< "\n\nmu:\n" << mu
			<< "\n\nsigma:\n" << sigma << endl;*/

		return NO_ERRORS;
	} // createNormalizedDataset

} // dmt

#endif // DATA_MINING_H


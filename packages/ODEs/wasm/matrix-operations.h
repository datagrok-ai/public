#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#include "../../../../Eigen/Eigen/Dense"
using namespace Eigen;

namespace matrOper {

	/* computation result code */
	enum ResultCode {
		NO_ERRORS = 0,
		UNKNOWN_PROBLEM		
	};

	/* Compute inverse matrix.
	     sourcePtr - source square matrix pointer
		 rowCount - number of rows
		 inversePtr - inverse matrix pointer */
	template<typename Type>
	int computeInverseMatrix(Type* sourcePtr, const int rowCount, Type* inversePtr) noexcept {
		Map<Matrix<Type, Dynamic, Dynamic, RowMajor>> A(sourcePtr, rowCount, rowCount);
		Map<Matrix<Type, Dynamic, Dynamic, RowMajor>> invA(inversePtr, rowCount, rowCount);

		invA = A.inverse();

		return NO_ERRORS;
	}
} // matrOper

#endif // !MATRIX_OPERATIONS_H


#include "matrix-operations.h"

#include <emscripten.h>

// The following provides convenient naming of the exported functions.
extern "C" {
	int invMatrixD(double* sourcePtr, const int rowCount, double* inversePtr) noexcept;
	int invMatrixF(float* sourcePtr, const int rowCount, float* inversePtr) noexcept;
};

EMSCRIPTEN_KEEPALIVE
/* Compute inverse matrix: double type. */
int invMatrixD(double* sourcePtr, const int rowCount, double* inversePtr) noexcept {
	return matrOper::computeInverseMatrix(sourcePtr, rowCount, inversePtr);
}

EMSCRIPTEN_KEEPALIVE
/* Compute inverse matrix: float type. */
int invMatrixF(float* sourcePtr, const int rowCount, float* inversePtr) noexcept {
	return matrOper::computeInverseMatrix(sourcePtr, rowCount, inversePtr);
}
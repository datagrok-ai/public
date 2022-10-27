#include <emscripten.h>

#include "..\..\Eigen\Eigen\Dense"

using namespace Eigen;

extern "C" {
    void sumByEigen(int * a, int aLength, int * b, int bLength, int * sum, int sumLength);
}

//name: sumOfColumns (non-export)
//input: column col1
//input: column col2
//output: column result [sum new(col1.length)]
EMSCRIPTEN_KEEPALIVE
void sumByEigen(int * a, int aLength, int * b, int bLength, int * sum, int sumLength)
{
	Map<RowVectorXi> vA(a, aLength);
	Map<RowVectorXi> vB(b, bLength);
	Map<RowVectorXi> vC(sum, sumLength);
	vC = vA + vB;
}
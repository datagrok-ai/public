// Test.cpp

#include <cmath>
using namespace std;

#include <emscripten.h>

extern "C" {
     int solveTest(float initial, float final, float step,
        float _xInitial, float _yInitial, 
        
        int _tCount, int _varsCount,
        float * _solution, int _solutionRowCount, int _solutionColCount) noexcept;
}

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "odes.h"
using namespace odes;

namespace Test
{
    // tollerance
    float _tol = 0.00005f;

    // dimension of solution
    const int DIM = 2;

    // parameters

    //the right part of the ODEs
    VectorXd _f(double t, VectorXd & _y) noexcept
    {
        VectorXd _res(DIM);

        // extract variables
        double x = _y(0);
        double y = _y(1);

        // output computation
        _res(0) =  cos(t);
        _res(1) =  - sin(t);

        return _res;
    } // _f

    // Jacobian (it is required, when applying implicit method)
    MatrixXd _J(double t, VectorXd & _y, double _eps) noexcept
    {
        MatrixXd _res(DIM, DIM);
        VectorXd _val = _f(t, _y);
        VectorXd _yDer = _y;

        for (int i = 0; i < DIM; i++) {
            _yDer(i) += _eps;
            _res.col(i) = (_f(t, _yDer) - _val) / _eps;
            _yDer(i) -= _eps;
        }

        return _res;
    } // _J

    // Derivative with respect to t (it is required, when applying implicit method)
    VectorXd _T(double t, VectorXd & _y, double _eps) noexcept
    {
        return (_f(t + _eps, _y) - _f(t, _y)) / _eps;
    } // _T

}; // Test

//name: solveTest
//input: double initial = 0.0 {caption: initial; category: time, minutes}
//input: double final = 5.0 {caption: final; category: time, minutes}
//input: double step = 0.1 {caption: step; category: time, minutes}
//input: double _xInitial = 0.0 {units: fhfhf; caption: x; category: initial values}
//input: double _yInitial = 1.0 {units: adasda; caption: y; category: initial values}
//input: int _tCount
//input: int _varsCount
//output: column_list _solution [new(_tCount, _varsCount); 't'; 'x(t)'; 'y(t)']
//output: dataframe dfSolution [_solution] {caption: Solution; viewer: Line chart(x: "t, time (minutes)", sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter") | Grid(block: 100) }
//editor: Compute:RichFunctionViewEditor
EMSCRIPTEN_KEEPALIVE
int solveTest(float initial, float final, float step,
    float _xInitial, float _yInitial, 
    
    int _tCount, int _varsCount,
    float * _solution, int _solutionRowCount, int _solutionColCount) noexcept
{
    using namespace Test;

    float _initialVals[DIM];

    // initial values
    _initialVals[0] = _xInitial;
    _initialVals[1] = _yInitial;

    // parameters

    return solveODE(_f, initial, final, step, _initialVals, _tol,
            _solution, _tCount, _varsCount);
} //solveTest
// TestExample.cpp

#include <cmath>
using namespace std;

#include <emscripten.h>

extern "C" {
     int solveTestExample(float initial, float final, float step,
        float _xInitial, float _yInitial, 
        float _param1Val, float _param2Val, 
        int _tCount, int _varsCount,
        float * _solution, int _solutionRowCount, int _solutionColCount) noexcept;
}

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "odes.h"
using namespace odes;

namespace TestExample
{
    // tollerance
    float _tol = 0.00005f;

    // dimension of solution
    const int DIM = 2;

    // parameters
    double param1 = 0.0;
    double param2 = 0.0;

    //the right part of the ODEs
    VectorXd _f(double t, VectorXd & _y) noexcept
    {
        VectorXd _res(DIM);

        // constants
        const double const1 = 1.0;
        const double const2 = 3.0;

        // extract variables
        double x = _y(0);
        double y = _y(1);

        // expressions
        double coef1 = const1 + param1;
        double coef2 = const2 + param2 + 0.0;

        // output computation
        _res(0) = coef1 * y;
        _res(1) = coef2 * x;

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

}; // TestExample

//name: solveTestExample
//input: double initial = 0.0 {caption: initial; category: time, minutes}
//input: double final = 5.0 {caption: final; category: time, minutes}
//input: double step = 0.1 {caption: step; category: time, minutes}
//input: double _xInitial = 2.0 {units: x y.o.; caption: x; category: initial values}
//input: double _yInitial = 0.0 {units: y y.o.; caption: y; category: initial values}
//input: double _param1Val = 1.0 {units: param1 y.o.; caption: param1; category: parameters}
//input: double _param2Val = -1.0 {units: param2 y.o.; caption: param2; category: parameters}
//input: int _tCount
//input: int _varsCount
//output: column_list _solution [new(_tCount, _varsCount); 't'; 'x(t)'; 'y(t)']
//output: dataframe dfSolution [_solution] {caption: Solution; viewer: Line chart(x: "t, time (minutes)", sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter") | Grid(block: 100) }
//editor: Compute:RichFunctionViewEditor
EMSCRIPTEN_KEEPALIVE
int solveTestExample(float initial, float final, float step,
    float _xInitial, float _yInitial, 
    float _param1Val, float _param2Val, 
    int _tCount, int _varsCount,
    float * _solution, int _solutionRowCount, int _solutionColCount) noexcept
{
    using namespace TestExample;

    float _initialVals[DIM];

    // initial values
    _initialVals[0] = _xInitial;
    _initialVals[1] = _yInitial;

    // parameters
    param1 = _param1Val;
    param2 = _param2Val;

    return solveODE(_f, _T, _J, initial, final, step, _initialVals, _tol,
            _solution, _tCount, _varsCount);
} //solveTestExample
// OneCompartmentPK.cpp

#include <cmath>
using namespace std;

#include <emscripten.h>

extern "C" {
     int solveOneCompartmentPK(float initial, float final, float step,
        float _depotInitial, float _centrInitial, float _periInitial, float _effInitial, 
        float _KAVal, float _CLVal, float _V2Val, float _QVal, float _V3Val, float _KinVal, float _KoutVal, float _EC50Val, 
        int _tCount, int _varsCount,
        float * _solution, int _solutionRowCount, int _solutionColCount) noexcept;
}

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "odes.h"
using namespace odes;

namespace OneCompartmentPK
{
    // tollerance
    float _tol = 0.00005f;

    // dimension of solution
    const int DIM = 4;

    // parameters
    double KA = 0.0;
    double CL = 0.0;
    double V2 = 0.0;
    double Q = 0.0;
    double V3 = 0.0;
    double Kin = 0.0;
    double Kout = 0.0;
    double EC50 = 0.0;

    //the right part of the ODEs
    VectorXd _f(double t, VectorXd & _y) noexcept
    {
        VectorXd _res(DIM);

        // extract variables
        double depot = _y(0);
        double centr = _y(1);
        double peri = _y(2);
        double eff = _y(3);

        // expressions
        double C2 =  centr / V2;
        double C3 =  peri / V3;

        // output computation
        _res(0) =  - KA * depot;
        _res(1) =  KA * depot - CL * C2;
        _res(2) =  Q * C2 - Q * C3;
        _res(3) =  Kin - Kout * (1.0 - C2 / (EC50 + C2)) * eff;

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

}; // OneCompartmentPK

//name: solveOneCompartmentPK
//input: double initial = 0.0 {caption: initial; category: time, hours}
//input: double final = 12.0 {caption: final; category: time, hours}
//input: double step = 0.01 {caption: step; category: time, hours}
//input: double _depotInitial = 0.0 {units: ; caption: depot; category: initial values}
//input: double _centrInitial = 0.0 {units: ; caption: centr; category: initial values}
//input: double _periInitial = 0.0 {units: ; caption: peri; category: initial values}
//input: double _effInitial = 1.0 {units: ; caption: eff; category: initial values}
//input: double _KAVal = 0.3 {units: ; caption: KA; category: parameters}
//input: double _CLVal = 2.0 {units: ; caption: CL; category: parameters}
//input: double _V2Val = 4.0 {units: ; caption: V2; category: parameters}
//input: double _QVal = 1.0 {units: ; caption: Q; category: parameters}
//input: double _V3Val = 30.0 {units: ; caption: V3; category: parameters}
//input: double _KinVal = 0.2 {units: ; caption: Kin; category: parameters}
//input: double _KoutVal = 0.2 {units: ; caption: Kout; category: parameters}
//input: double _EC50Val = 8.0 {units: ; caption: EC50; category: parameters}
//input: int _tCount
//input: int _varsCount
//output: column_list _solution [new(_tCount, _varsCount); 't'; 'depot(t)'; 'centr(t)'; 'peri(t)'; 'eff(t)']
//output: dataframe dfSolution [_solution] {caption: Solution; viewer: Line chart(x: "t", sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter") | Grid(block: 100) }
//editor: Compute:RichFunctionViewEditor
EMSCRIPTEN_KEEPALIVE
int solveOneCompartmentPK(float initial, float final, float step,
    float _depotInitial, float _centrInitial, float _periInitial, float _effInitial, 
    float _KAVal, float _CLVal, float _V2Val, float _QVal, float _V3Val, float _KinVal, float _KoutVal, float _EC50Val, 
    int _tCount, int _varsCount,
    float * _solution, int _solutionRowCount, int _solutionColCount) noexcept
{
    using namespace OneCompartmentPK;

    float _initialVals[DIM];

    // initial values
    _initialVals[0] = _depotInitial;
    _initialVals[1] = _centrInitial;
    _initialVals[2] = _periInitial;
    _initialVals[3] = _effInitial;

    // parameters
    KA = _KAVal;
    CL = _CLVal;
    V2 = _V2Val;
    Q = _QVal;
    V3 = _V3Val;
    Kin = _KinVal;
    Kout = _KoutVal;
    EC50 = _EC50Val;

    return solveODE(_f, _T, _J, initial, final, step, _initialVals, _tol,
            _solution, _tCount, _varsCount);
} //solveOneCompartmentPK

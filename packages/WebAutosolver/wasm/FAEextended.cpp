// FAEextended.cpp

#include <cmath>
using namespace std;

#include <emscripten.h>

extern "C" {
     int solveFAEextended(float initial, float final, float step,
        float _FFoxInitial, float _KKoxInitial, float _FFredInitial, float _KKredInitial, float _FfreeInitial, float _KfreeInitial, float _FKredInitial, float _FKoxInitial, float _MEAthiolInitial, float _CO2Initial, float _yO2PInitial, float _CystamineInitial, float _VLInitial, 
        float _qinVal, float _percentO2saturationVal, float _yO2inVal, float _pKa2MEAVal, float _HVal, float _TVal, float _RVal, float _PVal, float _TimeToSwitchVal, 
        int _tCount, int _varsCount,
        float * _solution, int _solutionRowCount, int _solutionColCount) noexcept;
}

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "odes.h"
using namespace odes;

namespace FAEextended
{
    // tollerance
    float _tol = 0.00005f;

    // dimension of solution
    const int DIM = 13;

    // parameters
    double qin = 0.0;
    double percentO2saturation = 0.0;
    double yO2in = 0.0;
    double pKa2MEA = 0.0;
    double H = 0.0;
    double T = 0.0;
    double R = 0.0;
    double P = 0.0;
    double TimeToSwitch = 0.0;

    //the right part of the ODEs
    VectorXd _f(double t, VectorXd & _y) noexcept
    {
        VectorXd _res(DIM);

        // constants
        const double VLinitial = 6.6;
        const double Vtotalvessel = 10.0;
        const double AgitatorSpeed = 400.0;
        const double AgitatorDiameter = 6.0;
        const double AgitatorPowerNumber = 2.1;
        const double pH = 7.4;
        const double k1red = 0.05604;
        const double k1ox = 0.0108;
        const double k2Fd = 1.35;
        const double k2Fa = 110400000.0;
        const double k2Kd = 0.04038;
        const double k2Ka = 120000000.0;
        const double k3FKa = 181200000.0;
        const double k3FKd = 0.01188;
        const double k4ox = 0.0108;
        const double k4red = 0.05604;
        const double kthiolox = 0.04;
        const double krcyst = 0.0;

        // extract variables
        double FFox = _y(0);
        double KKox = _y(1);
        double FFred = _y(2);
        double KKred = _y(3);
        double Ffree = _y(4);
        double Kfree = _y(5);
        double FKred = _y(6);
        double FKox = _y(7);
        double MEAthiol = _y(8);
        double CO2 = _y(9);
        double yO2P = _y(10);
        double Cystamine = _y(11);
        double VL = _y(12);

        // expressions
        double constForklasurface =  (3.932 * pow((pow(AgitatorSpeed, 3.0) * pow(AgitatorDiameter, 5.0) * AgitatorPowerNumber / 2160000000000.0), 0.361)) / 60.0;
        double klasurface =  pow(VL, - 0.65) * constForklasurface;
        double MEAthiolate =  MEAthiol * pow(10.0,(pH - pKa2MEA));
        double qout =  qin - klasurface * (yO2P * H - CO2) * VL * R * T / (P * 1000.0);
        double OTR =  klasurface * (yO2P * H - CO2);
        double Vg =  Vtotalvessel - VL;
        double Fin =  t < TimeToSwitch ? (0.0) : (0.025);
        double Fpermeate =  t < TimeToSwitch ? (0.025) : (Fin);
        double CO2in =  percentO2saturation * 7.17 / (32.0 * 100.0);
        double Vres =  VLinitial / VL;
        double MEAthiolate_t_by_Vres_t_squared =  pow(MEAthiolate * Vres, 2.0);
        double FFox_to_FFred =  k1red * FFox * Vres * MEAthiolate_t_by_Vres_t_squared;
        double FFred_to_FFox =  k1ox * FFred * Vres;
        double FFred_to_Ffree =  k2Fd * FFred * Vres;
        double Ffree_to_FFred =  k2Fa * pow(Ffree * Vres, 2.0) * MEAthiolate_t_by_Vres_t_squared;
        double KKox_to_KKred =  k1red * KKox * Vres * MEAthiolate_t_by_Vres_t_squared;
        double KKred_to_KKox =  k1ox * KKred * Vres;
        double KKred_to_Kfree =  k2Kd * KKred * Vres;
        double Kfree_to_KKred =  k2Ka * pow(Kfree * Vres, 2.0) * MEAthiolate_t_by_Vres_t_squared;
        double free_to_FKred =  k3FKa * Ffree * Vres * Kfree * Vres;
        double FKred_to_free =  k3FKd * FKred * Vres;
        double FKred_to_FKox =  k4ox * FKred * Vres * pow(Cystamine * Vres, 2.0);
        double FKox_to_FKred =  k4red * FKox * Vres * MEAthiolate_t_by_Vres_t_squared;
        double Vres_t_by_CO2 =  Vres * CO2;
        double sqrt_of_Vres_t_by_CO2 =  (Vres_t_by_CO2 >= 0.0) ? sqrt(Vres_t_by_CO2) : 0.0;
        double MEAthiol_t_by_Vres_t_squared =  pow(MEAthiol * Vres, 2.0);

        // output computation
        _res(0) =  - FFox_to_FFred + FFred_to_FFox;
        _res(1) =  - KKox_to_KKred + KKred_to_KKox;
        _res(2) =  FFox_to_FFred - FFred_to_FFox - FFred_to_Ffree + Ffree_to_FFred;
        _res(3) =  KKox_to_KKred - KKred_to_KKox - KKred_to_Kfree + Kfree_to_KKred;
        _res(4) =  2.0 * FFred_to_Ffree - 2.0 * Ffree_to_FFred - free_to_FKred + FKred_to_free;
        _res(5) =  2.0 * KKred_to_Kfree - 2.0 * Kfree_to_KKred - free_to_FKred + FKred_to_free;
        _res(6) =  free_to_FKred - FKred_to_free - FKred_to_FKox + FKox_to_FKred;
        _res(7) =  FKred_to_FKox - FKox_to_FKred;
        _res(8) =  2.0 * ( - FFox_to_FFred + FFred_to_FFox - KKox_to_KKred + KKred_to_KKox + FFred_to_Ffree + KKred_to_Kfree - Ffree_to_FFred - Kfree_to_KKred - FKox_to_FKred - kthiolox * MEAthiol_t_by_Vres_t_squared * sqrt_of_Vres_t_by_CO2) - (MEAthiol + MEAthiolate) * (Fin + Fpermeate) / VL;
        _res(9) =  (Fin * CO2in - 2.0 * Fpermeate * CO2) / VL + OTR - 0.5 * kthiolox * MEAthiol_t_by_Vres_t_squared * sqrt_of_Vres_t_by_CO2;
        _res(10) =  - OTR * (VL / Vg) * R * T * P + yO2in * qin - yO2P * qout;
        _res(11) =  kthiolox * MEAthiol_t_by_Vres_t_squared * sqrt_of_Vres_t_by_CO2 - krcyst * Cystamine * Vres - (Fin + Fpermeate) * Cystamine / VL;
        _res(12) =  Fin - Fpermeate;

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

}; // FAEextended

//name: solveFAEextended
//input: double initial = 0.0 {caption: initial; category: time, minutes}
//input: double final = 1000.0 {caption: final; category: time, minutes}
//input: double step = 0.1 {caption: step; category: time, minutes}
//input: double _FFoxInitial = 0.268 {units: ; caption: FFox; category: initial values}
//input: double _KKoxInitial = 0.268 {units: ; caption: KKox; category: initial values}
//input: double _FFredInitial = 0.0 {units: ; caption: FFred; category: initial values}
//input: double _KKredInitial = 0.0 {units: ; caption: KKred; category: initial values}
//input: double _FfreeInitial = 0.0 {units: ; caption: Ffree; category: initial values}
//input: double _KfreeInitial = 0.0 {units: ; caption: Kfree; category: initial values}
//input: double _FKredInitial = 0.0 {units: ; caption: FKred; category: initial values}
//input: double _FKoxInitial = 0.0 {units: ; caption: FKox; category: initial values}
//input: double _MEAthiolInitial = 34.0 {units: ; caption: MEAthiol; category: initial values}
//input: double _CO2Initial = 0.22 {units: ; caption: CO2; category: initial values}
//input: double _yO2PInitial = 0.209 {units: ; caption: yO2P; category: initial values}
//input: double _CystamineInitial = 0.0 {units: ; caption: Cystamine; category: initial values}
//input: double _VLInitial = 6.6 {units: ; caption: VL; category: initial values}
//input: double _qinVal = 1.0 {units: Liters gas/minute; caption: qin; category: parameters}
//input: double _percentO2saturationVal = 100.0 {units: ; caption: percentO2saturation; category: parameters}
//input: double _yO2inVal = 0.209 {units: ; caption: yO2in; category: parameters}
//input: double _pKa2MEAVal = 8.19 {units: ; caption: pKa2MEA; category: parameters}
//input: double _HVal = 1.072069378 {units: ; caption: H; category: parameters}
//input: double _TVal = 300.0 {units: degK; caption: T; category: parameters}
//input: double _RVal = 0.082 {units: Liter Atm / mole degK; caption: R; category: parameters}
//input: double _PVal = 1.0 {units: atma; caption: P; category: parameters}
//input: double _TimeToSwitchVal = 180.0 {units: ; caption: TimeToSwitch; category: parameters}
//input: int _tCount
//input: int _varsCount
//output: column_list _solution [new(_tCount, _varsCount); 't'; 'FFox(t)'; 'KKox(t)'; 'FFred(t)'; 'KKred(t)'; 'Ffree(t)'; 'Kfree(t)'; 'FKred(t)'; 'FKox(t)'; 'MEAthiol(t)'; 'CO2(t)'; 'yO2P(t)'; 'Cystamine(t)'; 'VL(t)']
//output: dataframe dfSolution [_solution] {caption: Solution; viewer: Line chart(x: "t, time (minutes)", sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter") | Grid(block: 100) }
//editor: Compute:RichFunctionViewEditor
EMSCRIPTEN_KEEPALIVE
int solveFAEextended(float initial, float final, float step,
    float _FFoxInitial, float _KKoxInitial, float _FFredInitial, float _KKredInitial, float _FfreeInitial, float _KfreeInitial, float _FKredInitial, float _FKoxInitial, float _MEAthiolInitial, float _CO2Initial, float _yO2PInitial, float _CystamineInitial, float _VLInitial, 
    float _qinVal, float _percentO2saturationVal, float _yO2inVal, float _pKa2MEAVal, float _HVal, float _TVal, float _RVal, float _PVal, float _TimeToSwitchVal, 
    int _tCount, int _varsCount,
    float * _solution, int _solutionRowCount, int _solutionColCount) noexcept
{
    using namespace FAEextended;

    float _initialVals[DIM];

    // initial values
    _initialVals[0] = _FFoxInitial;
    _initialVals[1] = _KKoxInitial;
    _initialVals[2] = _FFredInitial;
    _initialVals[3] = _KKredInitial;
    _initialVals[4] = _FfreeInitial;
    _initialVals[5] = _KfreeInitial;
    _initialVals[6] = _FKredInitial;
    _initialVals[7] = _FKoxInitial;
    _initialVals[8] = _MEAthiolInitial;
    _initialVals[9] = _CO2Initial;
    _initialVals[10] = _yO2PInitial;
    _initialVals[11] = _CystamineInitial;
    _initialVals[12] = _VLInitial;

    // parameters
    qin = _qinVal;
    percentO2saturation = _percentO2saturationVal;
    yO2in = _yO2inVal;
    pKa2MEA = _pKa2MEAVal;
    H = _HVal;
    T = _TVal;
    R = _RVal;
    P = _PVal;
    TimeToSwitch = _TimeToSwitchVal;

    return solveODE(_f, _T, _J, initial, final, step, _initialVals, _tol,
            _solution, _tCount, _varsCount);
} //solveFAEextended
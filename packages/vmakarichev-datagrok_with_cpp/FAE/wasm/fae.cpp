// fae.cpp

#include<cmath>
using std::pow;
using std::sqrt;
using std::abs;

#include <emscripten.h>

extern "C" {
    int solveFAE(float t0, float t1, float h, 
                 float FFox, float KKox, float FFred, float KKred,
                 float Ffree, float Kfree, float FKred, float FKox,
                 float MEAthiol_t, float CO2, float yO2P, float Cystamine,
                 float VL, int timesCount, int varsCount,
                 float * result, int resultRowCount, int resultColCount) noexcept;
}

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "odeSolver.h"
using namespace ode;

namespace fae
{
	// overall tollerance
    float tol = 0.00005f;

    // dimension of solution
	const int DIM = 13; 
   
	// Version 3
	VectorXd f(double _time, VectorXd & y) noexcept
	{
		VectorXd res(y.size());

		// ORIGINAL EXPRESSIONS ARE REMOVED!

		return res;
	} // f
 

}; //fae

//name: solveFAE
//input: double t0 = 0 {units: minutes; caption: initial; category: Time}
//input: double t1 = 1000 {units: minutes; caption: final; category: Time}
//input: double h = 0.01 {units: minutes; caption: step; category: Time}
//input: double FFox = 0.268 {caption: FFox; category: Initial values}
//input: double KKox = 0.268 {caption: KKox; category: Initial values}
//input: double FFred = 0.0 {caption: FFred; category: Initial values}
//input: double KKred = 0.0 {caption: KKred; category: Initial values}
//input: double Ffree = 0.0 {caption: Ffree; category: Initial values}
//input: double Kfree = 0.0 {caption: Kfree; category: Initial values}
//input: double FKred = 0.0 {caption: FKred; category: Initial values}
//input: double FKox = 0.0 {caption: FKox; category: Initial values}
//input: double MEAthiol = 34.0 {caption: MEAthiol; category: Initial values}
//input: double CO2 = 0.22 {caption: CO2; category: Initial values}
//input: double yO2P = 0.209 {caption: yO2P; category: Initial values}
//input: double Cystamine = 0.0 {caption: Cystamine; category: Initial values}
//input: double VL = 6.6 {caption: VL; category: Initial values}
//input: int timesCount
//input: int varsCount
//output: column_list result [new(timesCount, varsCount); 't, time (minutes)'; 'FFox(t)'; 'KKox(t)'; 'FFred(t)'; 'KKred(t)'; 'Ffree(t)'; 'Kfree(t)'; 'FKred(t)'; 'FKox(t)'; 'MEAthiol(t)'; 'CO2(t)'; 'yO2P(t)'; 'Cystamine(t)'; 'VL(t)']
//output: dataframe solution [result] {caption: Solution; viewer: Line chart(x: "t, time (minutes)", sharex: "true", multiAxis: "true", yGlobalScale: "true", multiAxisLegendPosition: "RightCenter") | Grid(block: 100) }
//editor: Compute:RichFunctionViewEditor
EMSCRIPTEN_KEEPALIVE
int solveFAE(float t0, float t1, float h, 
  float FFox, float KKox, float FFred, float KKred,
  float Ffree, float Kfree, float FKred, float FKox,
  float MEAthiol, float CO2, float yO2P, float Cystamine,
  float VL, int timesCount, int varsCount,
  float * result, int resultRowCount, int resultColCount) noexcept
{
    using namespace fae;

	float yInitial[DIM];

	yInitial[0] = FFox;
	yInitial[1] = KKox;
	yInitial[2] = FFred;
	yInitial[3] = KKred;
	yInitial[4] = Ffree;
	yInitial[5] = Kfree;
	yInitial[6] = FKred;
	yInitial[7] = FKox;
	yInitial[8] = MEAthiol;
	yInitial[9] = CO2;
	yInitial[10] = yO2P;
	yInitial[11] = Cystamine;
	yInitial[12] = VL;

    return solveODE(f, t0, t1, h, yInitial, tol, result, timesCount, varsCount);
} // solveFAE


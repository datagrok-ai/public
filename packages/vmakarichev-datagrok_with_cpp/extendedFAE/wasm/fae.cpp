// fae.cpp

#include<cmath>
using std::pow;
using std::sqrt;
using std::abs;

#include <emscripten.h>

extern "C" {
    int fullFAE(float t0, float t1, float h, 
                float FFox, float KKox, float FFred, float KKred,
				float Ffree, float Kfree, float FKred, float FKox,
                float MEAthiol, float CO2, float yO2P, float Cystamine,
                float VL, float qinVal, float percentO2saturationVal,
                float yO2inVal, float pKa2MEAVal, float HVal, float TVal,
                float RVal, float PVal,
                int timesCount, int varsCount,
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

	// parameters
	double qin = 0.0;
	double percentO2saturation = 0.0;
	double yO2in = 0.0;
	double pKa2MEA = 0.0;
	double H = 0.0;
	double T = 0.0;
	double R = 0.0;
	double P = 0.0;    

	VectorXd f(double _time, VectorXd & y) noexcept
	{
		VectorXd res(y.size());		

		// ORIGINAL EXPRESSIONS ARE DELETED!
		
		return res;
	} // f
 

}; //fae

// ORIGINAL VALUES ARE DELETED!

//name: fullFAE
//input: double t0 = 0 {units: minutes; caption: initial; category: Time}
//input: double t1 =  {units: minutes; caption: final; category: Time}
//input: double h =  {units: minutes; caption: step; category: Time}
//input: double FFox =  {caption: FFox; category: Initial values}
//input: double KKox =  {caption: KKox; category: Initial values}
//input: double FFred =  {caption: FFred; category: Initial values}
//input: double KKred =  {caption: KKred; category: Initial values}
//input: double Ffree =  {caption: Ffree; category: Initial values}
//input: double Kfree =  {caption: Kfree; category: Initial values}
//input: double FKred =  {caption: FKred; category: Initial values}
//input: double FKox =  {caption: FKox; category: Initial values}
//input: double MEAthiol =  {caption: MEAthiol; category: Initial values}
//input: double CO2 =  {caption: CO2; category: Initial values}
//input: double yO2P =  {caption: yO2P; category: Initial values}
//input: double Cystamine =  {caption: Cystamine; category: Initial values}
//input: double VL =  {caption: VL; category: Initial values}
//input: double qinVal =  {caption: qin; category: Parameters}
//input: double percentO2saturationVal =  {caption: percentO2saturation; category: Parameters}
//input: double yO2inVal =  {caption: yO2in; category: Parameters}
//input: double pKa2MEAVal =  {caption: pKa2MEA; category: Parameters}
//input: double HVal =  {caption: H; category: Parameters}
//input: double TVal =  {caption: T; category: Parameters}
//input: double RVal =  {caption: R; category: Parameters}
//input: double PVal =  {caption: P; category: Parameters} 
//input: int timesCount
//input: int varsCount
//output: column_list result [new(timesCount, varsCount); 't, time (minutes)'; 'FFox(t)'; 'KKox(t)'; 'FFred(t)'; 'KKred(t)'; 'Ffree(t)'; 'Kfree(t)'; 'FKred(t)'; 'FKox(t)'; 'MEAthiol(t)'; 'CO2(t)'; 'yO2P(t)'; 'Cystamine(t)'; 'VL(t)']
//output: dataframe solution [result] {caption: Solution; viewer: Line chart(x: "t, time (minutes)", sharex: "true", multiAxis: "true", yGlobalScale: "true", multiAxisLegendPosition: "RightCenter") | Grid(block: 100) }
//editor: Compute:RichFunctionViewEditor
EMSCRIPTEN_KEEPALIVE
int fullFAE(float t0, float t1, float h, 
  float FFox, float KKox, float FFred, float KKred,
  float Ffree, float Kfree, float FKred, float FKox,
  float MEAthiol, float CO2, float yO2P, float Cystamine,
  float VL, float qinVal, float percentO2saturationVal,
  float yO2inVal, float pKa2MEAVal, float HVal, float TVal,
  float RVal, float PVal,
  int timesCount, int varsCount,
  float * result, int resultRowCount, int resultColCount) noexcept
{
    using namespace fae;

	float yInitial[DIM];

	// initial values
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

	// parameters
	qin = qinVal;
	percentO2saturation = percentO2saturationVal;
	yO2in = yO2inVal;
	pKa2MEA = pKa2MEAVal;
	H = HVal;
	T = TVal;
	R = RVal;
	P = PVal; 

    return solveODElight(f, t0, t1, h, yInitial, tol, result, timesCount, varsCount);
} // solveFAE


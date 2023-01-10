// fae.cpp

#include<cmath>
using std::pow;
using std::sqrt;

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
//input: double t0
//input: double t1
//input: double h
//input: double FFox
//input: double KKox
//input: double FFred
//input: double KKred
//input: double Ffree
//input: double Kfree
//input: double FKred
//input: double FKox
//input: double MEAthiol
//input: double CO2
//input: double yO2P
//input: double Cystamine
//input: double VL
//input: int timesCount
//input: int varsCount
//output: column_list result [new(timesCount, varsCount)]
//output: dataframe solution [result]
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


// odeSolver.h

// ODE solvers

/* The formulas applied are taken from

   [1] Endre Suli and David F. Mayers. An Introduction to Numerical Analysis, 2003.

   [2] Steve Chapra and Raymond P. Canale. Numerical Methods for Engineers, 2021.

   [3] https://doi.org/10.1137/S1064827594276424

   [4] https://doi.org/10.1016/S0898-1221(00)00175-9
*/

#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include <cstring>

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

namespace ode {

	enum ResultCode {NO_ERRORS = 0, UNKNOWN_PROBLEM };	

	/* One step of the Runge-Kutta method of the order 3 for the problem dy/dt = f(t,y), y(t0) = y0.
	   Returns value of y(t + h) computed using fourth-order R.-K. method.
	     f - right part of the differential equation;
	     t - current value of the independent variable;
	     y - solution at the point t;
	     h - step. */
	template<typename ArgType, typename VecType>
	VecType getNextPointRK3(VecType(*f)(ArgType, VecType &), ArgType t, VecType & y, ArgType h)
	{
		// The classic 3rd-order Runge-Kutta method
		VecType k1 = f(t, y);
		VecType yDer = y + 0.5 * h * k1;
		VecType k2 = f(t + h * 0.5, yDer);
		yDer = y - h * k1 + 2.0 * h * k2;
		VecType k3 = f(t + h * 0.5, yDer);		

		return y + (k1 + k2 * 4.0 + k3 ) * h / 6.0;
	}

	/* One step of the Runge-Kutta method of the order 4 for the problem dy/dt = f(t,y), y(t0) = y0.
	   Returns value of y(t + h) computed using fourth-order R.-K. method.
	      f - right part of the differential equation;
		  t - current value of the independent variable;
		  y - solution at the point t;
		  h - step. */
	template<typename ArgType, typename VecType>
	VecType getNextPointRK4(VecType (*f)(ArgType, VecType &), ArgType t, VecType & y, ArgType h)
	{
		// The classic fourth-order Runge-Kutta method
		VecType k1 = f(t, y);
		VecType yDer = y + 0.5 * h * k1;
		VecType k2 = f(t + h * 0.5, yDer);
		yDer = y + 0.5 * h * k2;
		VecType k3 = f(t + h * 0.5, yDer);
		yDer = y + h * k3;
		VecType k4 = f(t + h * 0.5, yDer);

		return y + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * h / 6.0;
	}

	/* One step of the Runge-Kutta method of the order 5 (Butcher's method) for the problem dy/dt = f(t,y), y(t0) = y0.
	   Returns value of y(t + h) computed using fourth-order R.-K. method.
	      f - right part of the differential equation;
		  t - current value of the independent variable;
		  y - solution at the point t;
		  h - step. */
	template<typename ArgType, typename VecType>
	VecType getNextPointRK5(VecType (*f)(ArgType, VecType &), ArgType t, VecType & y, ArgType h)
	{
		// Butcher’s (1964) fifth-order RK method
		VecType k1 = f(t, y);
		VecType yDer = y + 0.25 * h * k1;
		VecType k2 = f(t + h * 0.25, yDer);
		yDer = y + 0.125 * h * k1 + 0.125 * h * k2;
		VecType k3 = f(t + h * 0.25, yDer);
		yDer = y - 0.5 * h * k2 + h * k3;
		VecType k4 = f(t + h * 0.5, yDer);
		yDer = y + 0.1875 * h * k1 - 0.5625 * h * k4;
		VecType k5 = f(t + h * 0.75, yDer);
		yDer = y + (-3.0 * h * k1 + 2.0 * h * k2 + 12.0 * h * (k3 - k4) + 8.0 * h * k5 ) / 7.0;
		VecType k6 = f(t + h, yDer);

		return y + (7.0 * (k1 + k6) + 32.0 * (k3 + k5)  + 12.0 * k4) * h / 90.0;
	}

	/*  One-step ODE solver of the problem dy/dt = f(t, y), y( t0 ) = y0: one-dimensional case.
	    Returns NO_ERRORS in the case of success computations.
	      f - right part of the ODE solved;
		  times - array of times;
		  timesCount - length of the array times;
		  yInitial - initial value of the function y, i.e. y0 = y( times[0] );
		  getNextPoint - one-step method that provides value of the solution at the next point: it computes y(t + h) for given y(t);
		  solution - array for solution: solution[i] contains y( times[i] ).
	*/
	template<typename ArgType>
	int oneStepSolver(ArgType (*f)(ArgType, ArgType &), 
		ArgType *times, unsigned timesCount,
		ArgType yInitial,
		ArgType (*getNextPoint)(ArgType (*)(ArgType, ArgType &), ArgType, ArgType &, ArgType),
		ArgType *solution)
	{
		// solution at the initial point
		solution[0] = yInitial;

		// solution at the rest of points
		for (unsigned i = 1; i < timesCount; i++)
			solution[i] = getNextPoint(f, times[i - 1], solution[i - 1], times[i] - times[i - 1]);

		return NO_ERRORS;
	}

	/*  One-step ODE solver of the problem dy/dt = f(t, y), y( t0 ) = y0: multi-dimensional case.
	    Returns NO_ERRORS in the case of success computations.
	      f - right part of the ODE solved;
		  times - array of times;
		  timesCount - length of the array times;
		  yInitial - initial value of the function y, i.e. y0 = y( times[0] );
		  getNextPoint - one-step method that provides value of the solution at the next point: it computes y(t + h) for given y(t);
		  solution - array that contains solution.

	    REMARK. Solution of the current multi-dimension problem is a vector function y = y(t). 
		        Note that y(t) is a column vector for each t. Currently, solution is a matrix Y. 
				Each column of Y contains solution at the specific point t,	i.e. i-th column of Y is y( times[i] ).
				For the purpose of the further use in DATAGROK, the array "solution" is a concatanation of the matrix Y rows.
	*/
	template<typename ArgType, typename VecType>
	int oneStepSolver(VecType(*f)(ArgType, VecType &),
		ArgType *times, unsigned timesCount,
		ArgType * yInitial, unsigned dimCount,
		VecType(*getNextPoint)(VecType(*)(ArgType, VecType &), ArgType, VecType &, ArgType),
		ArgType *solution)
	{
		// matrix that contains solution: i-th column contains solution (a vector) at the point times[i]
		Map<Matrix<ArgType, Dynamic, Dynamic, RowMajor>> Y(solution, dimCount, timesCount);
				
		// vector for storing solution at the previos time 
		VecType yPrev(dimCount);

		// solution at the initial point: copying bytes apporach
		std::memcpy(yPrev.data(), yInitial, sizeof(ArgType) * dimCount);

		// store solution at the initial point
		Y.col(0) = yPrev;
		
		// solution at the rest of points
		for (unsigned i = 1; i < timesCount; i++)
		{
			yPrev = getNextPoint(f, times[i - 1], yPrev, times[i] - times[i - 1]);
			Y.col(i) = yPrev;
		}

		return NO_ERRORS;
	}


	/*  One-step ODE implicit solver of the problem dy/dt = f(t, y), y( t0 ) = y0: multi-dimensional case.
		Returns NO_ERRORS in the case of success computations.
		  f - right part of the ODE solved;
		  T - vector of derivatives df/dt, computed approximately;
		  J - matrix of derivatives df/dy, computed approximately;
		  times - array of times;
		  timesCount - length of the array times;
		  yInitial - initial value of the function y, i.e. y0 = y( times[0] );		  
		  solution - array that contains solution.

		REMARK. Solution of the current multi-dimension problem is a vector function y = y(t).
				Note that y(t) is a column vector for each t. Currently, solution is a matrix Y.
				Each column of Y contains solution at the specific point t,	i.e. i-th column of Y is y( times[i] ).
				For the purpose of the further use in DATAGROK, the array "solution" is a concatanation of the matrix Y rows.
	*/
	template<typename ArgType, typename VecType, typename MatType>
	int oneStepSolver(VecType(*f)(ArgType, VecType&),
		VecType(*T)(ArgType, VecType&, ArgType),
		MatType(*J)(ArgType, VecType&, ArgType),
		ArgType* times, unsigned timesCount,
		ArgType* yInitial, unsigned dimCount,		
		ArgType* solution)
	{
		// Value of h / eps, where h - step of ODEs solver, eps - step for computing derivatives
		const ArgType H_TO_EPS_RATIO = 0.1; 

		// matrix that contains solution: i-th column contains solution (a vector) at the point times[i]
		Map<Matrix<ArgType, Dynamic, Dynamic, RowMajor>> Y(solution, dimCount, timesCount);

		// vector for storing solution at the current time 
		VecType y(dimCount);

		// solution at the initial point: copying bytes apporach
		std::memcpy(y.data(), yInitial, sizeof(ArgType) * dimCount);

		// store solution at the initial point
		Y.col(0) = y;

		// Quantities used in Rosenbrock method (see papers [3, 4] for more details)
		ArgType d = static_cast<ArgType>(1.0 - sqrt(2.0) / 2.0);
		MatType I = MatType::Identity(dimCount, dimCount);
		MatType W(dimCount, dimCount);
		MatType invW(dimCount, dimCount);
		VecType f0(dimCount);
		VecType k1(dimCount);
		VecType f1(dimCount);
		VecType k2(dimCount);
		VecType yDer(dimCount);		

		// solution at the rest of points
		for (unsigned i = 1; i < timesCount; i++)
		{
			ArgType h = times[i] - times[i - 1];
			ArgType eps = h * H_TO_EPS_RATIO;
			ArgType t = times[i - 1];
						
			f0 = h * f(t, y);			
			W = I - h * d * J(t, y, eps); 
			invW = W.inverse();
			k1 = invW * (f0 + h * d * T(t, y, eps));
			yDer = y + 0.5 * h * k1;
			f1 = f(t + 0.5 * h, yDer);
			k2 = invW * (f1 - k1) + k1;
			yDer = y;
			y += k2 * h;
			Y.col(i) = y;

		} // for

		return NO_ERRORS;
	}


}; // ode


#endif


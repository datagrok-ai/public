// odeSolver.h

// ODE solvers

/* The formulas applied are taken from 
   [1] Endre Suli and David F. Mayers. An Introduction to Numerical Analysis, Cambridge University Press, 2003.
   
*/

#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include <cstring>

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

namespace ode {

	enum ResultCode {NO_ERRORS = 0, UNKNOWN_PROBLEM };

	/* One step of the Runge-Kutta method of the order 4 for the problem dy/dt = f(t,y).
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

	template<typename ArgType>
	int oneStepExplicitSolver(ArgType (*f)(ArgType, ArgType &), 
		ArgType *times, unsigned timesCount,
		ArgType yInitial,
		ArgType (*getNextPoint)(ArgType (*)(ArgType, ArgType &), ArgType, ArgType &, ArgType),
		ArgType *solution)
	{
		solution[0] = yInitial;

		for (unsigned i = 1; i < timesCount; i++)
			solution[i] = getNextPoint(f, times[i - 1], solution[i - 1], times[i] - times[i - 1]);

		return NO_ERRORS;
	}

	template<typename ArgType, typename VecType>
	int oneStepExplicitSolver(VecType(*f)(ArgType, VecType &),
		ArgType *times, unsigned timesCount,
		ArgType * yInitial, unsigned dimCount,
		VecType(*getNextPoint)(VecType(*)(ArgType, VecType &), ArgType, VecType &, ArgType),
		ArgType *solution)
	{
		Map<Matrix<ArgType, Dynamic, Dynamic, ColMajor>> Y(solution, dimCount, timesCount);

		std::memcpy(solution, yInitial, sizeof(ArgType) * dimCount);
		
		VecType yPrev(dimCount);

		for (unsigned i = 1; i < timesCount; i++)
		{
			yPrev = Y.col(i - 1);
			Y.col(i) = getNextPoint(f, times[i - 1], yPrev, times[i] - times[i - 1]);
		}

		return NO_ERRORS;
	}


}; // ode


#endif


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
#include <vector>
//#include <list>
//#include <deque>

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

namespace ode {

	enum ResultCode {NO_ERRORS = 0, UNKNOWN_PROBLEM, METHOD_FAILS };	

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

	/*  Fixed-step ODE solver of the problem dy/dt = f(t, y), y( t0 ) = y0: one-dimensional case.
	    Returns NO_ERRORS in the case of success computations.
	      f - right part of the ODE solved;
		  times - array of times;
		  timesCount - length of the array times;
		  yInitial - initial value of the function y, i.e. y0 = y( times[0] );
		  getNextPoint - one-step method that provides value of the solution at the next point: it computes y(t + h) for given y(t);
		  solution - array for solution: solution[i] contains y( times[i] ).
	*/
	template<typename ArgType>
	int fixedStepSolver(ArgType (*f)(ArgType, ArgType &), 
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

	/*  Fixed-step ODE solver of the problem dy/dt = f(t, y), y( t0 ) = y0: multi-dimensional case.
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
	int fixedStepSolver(VecType(*f)(ArgType, VecType &),
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


	/*  Fixed-step ODE implicit solver of the problem dy/dt = f(t, y), y( t0 ) = y0: multi-dimensional case.
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
	int fixedStepSolver(VecType(*f)(ArgType, VecType&),
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

	/*  One stage of Runge-Kutta-Cash-Karp method for solving the problem dy/dt = f(t,y), y(t0) = y0. 
	    Computes yOut and yErr, which are respectively numerical solution and errors at the point t + h.
		  f - right part of the ODE solved;
		  y - solution at the point t;
		  dydt - derivative of the solution at the point t;
		  t - current value of the independent variable;
		  h - step;
		  yOut - numerical solution at the point t + h;
		  yErr - errors at the point t + h. 
	*/
	template<typename ArgType, typename VecType>
	int RKCK( VecType(*f)(ArgType, VecType&), 
		VecType& y, VecType& dydt, ArgType t, ArgType h, VecType& yOut, VecType& yErr )
	{
		// implementation of R.-K. Cash-Karp method (see [2] for more details)
		VecType ytemp = y + 0.2 * h * dydt;
		VecType k2 = f(t + 0.2 * h, ytemp);
		ytemp = y + h * (0.075 * dydt + 0.225 * k2);
		VecType k3 = f(t + 0.3 * h, ytemp);
		ytemp = y + h * (0.3 * dydt - 0.9 * k2 + 1.2 * k3);
		VecType k4 = f(t + 0.6 * h, ytemp);
		ytemp = y + h * (-0.2037037037037037 * dydt + 2.5 * k2 - 2.592592592592593 * k3 + 1.296296296296296 * k4);
		VecType k5 = f(t + 0.6 * h, ytemp);
		ytemp = y + h * (0.0294958043981481 * dydt + 0.341796875 * k2 + 0.0415943287037037 * k3
			+ 0.4003454137731481 * k4 + 0.061767578125 * k5);
		VecType k6 = f(t + 0.875 * h, ytemp);
		yOut = y + h * (0.0978835978835979 * dydt + 0.4025764895330113 * k3
			+ 0.21043771043771045 * k4 + 0.2891022021456804 * k6);
		yErr = h * (-0.004293774801587311 * dydt + 0.018668586093857853 * k3 - 0.034155026830808066 * k4
			- 0.019321986607142856 * k5 + 0.03910220214568039 * k6);

		return NO_ERRORS;
	} // RKCK

	/*  One stage of Runge-Kutta-Cash-Karp method for solving the problem dy/dt = f(t,y), y(t0) = y0.
		Computes an appropriate step h, which provides the desired error, 
		solution at the corresponding point and length of the next step.
		  f - right part of the ODE solved;
		  y - solution at the point t;
		  dydt - derivative of the solution at the point t;
		  t - current value of the independent variable;
		  hTry - initial value of the step;
		  yScale - determines how the error is scaled;
		  hNext - initial step for computing solution at the next point;
		  tol - overall tolerance level.
	*/
	template<typename ArgType, typename VecType>
	int adaptiveRKCK(VecType(*f)(ArgType, VecType&),
		VecType& y, VecType& dydt, ArgType& t, ArgType hTry, VecType& yScale, ArgType& hNext, ArgType tol )
	{
		// implementation of one adaptive step of R.-K. Cash-Karp method (see [2] for more details)

		// hyperparameters of the method
		const ArgType SAFETY = 0.9;
		const ArgType PSHRNK = -0.25;
		const ArgType PSGROW = -0.2;
		const ArgType REDUCE_COEF = 0.25;
		const ArgType GROW_COEF = 4.0;
		const ArgType ERR_CONTR = 1.89e-4;

		// initialization
		ArgType h = hTry;
		unsigned dim = y.size();
		VecType yTemp(dim);
		VecType yErr(dim);

		// computation of solution (y), time (t) and next step (hNext)
		while (true)
		{
			// apply one step of R.-K. Cash-Karp computations
			int resultCode = RKCK(f, y, dydt, t, h, yTemp, yErr);

			// check result code
			if (resultCode != NO_ERRORS)
				return resultCode;

			// estimating error
			ArgType errmax = 0.0;
			for (unsigned i = 0; i < dim; i++)
				errmax = fmax(errmax, fabs(yErr(i) / yScale(i)));
			errmax /= tol;

			// processing the error obtained
			if (errmax > static_cast<ArgType>(1)) {
				ArgType hTemp = SAFETY * h * pow(errmax, PSHRNK);
				h = fmax(hTemp, REDUCE_COEF * h);
				ArgType tNew = t + h;
				if (tNew == t)
					return METHOD_FAILS;
			} // if
			else {
				if (errmax > ERR_CONTR)
					hNext = SAFETY * h * pow(errmax, PSGROW);
				else
					hNext = GROW_COEF * h;
				t = t + h;
				y = yTemp;
				break;
			} // else
		} // while

		return NO_ERRORS;
	} // adaptiveRKCK

	/*  Explicit adaptive step ODEs solver for the problem dy/dt = f(t,y), y(t0) = y0.
	    Solves the initial problem on the segment [t0, t1] and stores results in the structures 
		times (independent variable t) and solutions (numerical solution at the correspondent t).
		  f - right part of the equation;
		  t0 - initial point;
		  t1 - end point;
		  hInitial - initial value of the step;
		  y0 - initial value of the solution, i.e. y(t0);
		  tol - overall tolerance level;
		  times - data structure for values of independent variable, i.e. t;
		  solutions - data structure for storing numerical solutions, i.e. y(t).
	*/
	template<typename ArgType, typename VecType, class ArgStruct, class VecStruct>
	int RKCKsolver(VecType(*f)(ArgType, VecType&),
		ArgType t0, ArgType t1, ArgType hInitial, VecType & y0, ArgType tol, 
		ArgStruct & times, VecStruct & solutions)
	{
		// Here, Cash-Karp method is implemented (see [2] for more details).
		
		// store initial values
		times.push_back(t0);
		solutions.push_back(y0);

		// dimension of the solution 
		unsigned dim = y0.size();

		// method routine
		VecType tiny = VecType::Ones(dim) * 1e-20;
		ArgType t = t0;
		ArgType h = hInitial;
		ArgType hNext = 0.0;
		VecType y = y0;
		VecType dydt(dim);
		VecType yScale(dim);
		bool flag = true;

		// compute numerical solution
		while (flag)
		{
			// compute derivative
			dydt = f(t, y);
			
			// compute scale vector
			yScale = y.cwiseAbs() + (dydt * h).cwiseAbs() + tiny;

			// check end point
			if (t + h > t1) {
				h = t1 - t;
				flag = false;
			} 

			// perform one step of R.-K. Cash-Karp method
			int resultCode = adaptiveRKCK(f, y, dydt, t, h, yScale, hNext, tol);

			// check result of one step
			if (resultCode != NO_ERRORS)
				return resultCode;

			// store results
			times.push_back(t);
			solutions.push_back(y);

			h = hNext;
		} // while

		return NO_ERRORS;
	} // adaptiveStepSolver

	/*  Linear interpolant.
	    Provides linear interpolation of ODE's solution stored in structures times & solutions, i.e. {t} and {y(t)}, 
		into dataframe. Values are obtained on the segment [t0, t1] with the step h.
		  t0 - initial point;
		  t1 - end point;
		  h - step;
		  dataFrame - the dataframe computed;
		  rowCount - number of rows;
		  colCount - number of columns;
		  times - structure that contains values of the independent variable t, i.e. {t};
		  solutions - structure that contains values of y(t).
	*/
	template<typename OperatingType, typename DataType, class ArgStruct, class VecStruct>
	int linearInterpolation(OperatingType t0, OperatingType t1, OperatingType h, 
		DataType * dataFrame, int rowCount, int colCount,
		ArgStruct & times, VecStruct & solutions)
	{
		// routine values
		auto tIter = times.begin();
		auto yIter = solutions.begin();

		auto tIterNext = times.begin();
		auto yIterNext = solutions.begin();
		++tIterNext;
		++yIterNext;

		OperatingType tLeft = 0.0;
		OperatingType tRight = 0.0;
		OperatingType cLeft = 0.0;
		OperatingType cRight = 0.0;
		OperatingType fLeft = 0.0;
		OperatingType fRight = 0.0;

		int index = 0;

		// the classic linear interpolation method
		for (OperatingType t = t0; t < t1; t += h)
		{
			while (t > *tIterNext) {
				++tIter;
				++yIter;
				++tIterNext;
				++yIterNext;
			}

			tLeft = *tIter;
			tRight = *tIterNext;

			cLeft = (tRight - t) / (tRight - tLeft);
			cRight = 1.0 - cLeft;

			dataFrame[index] = static_cast<DataType>(t);

			for (int j = 1; j < colCount; j++)
			{
				fLeft = (*yIter)(j - 1);
				fRight = (*yIterNext)(j - 1);

				dataFrame[index + j * rowCount] = static_cast<DataType>(cLeft * fLeft + cRight * fRight);
			}

			index++;
		}
		return NO_ERRORS;
	} // linearInterpolation

	/*  Solver of the initial ODE problem dy/dt = f(t,y), y(t0) = y0.
	    Solves the problem on the segment [t0, t1] and linearly interpolated with the step h results are stored in the dataframe.
		  f - right part of the equation;
		  t0 - initial point;
		  t1 - end point;
		  h - step;
		  y0 - initial value of the solution, i.e. y(t0);
		  tol - overall tolerance level;
		  dataframe - array that contains values of (t, y(t));
		  rowCount - number of rows;
		  colCount - number of columns.

		REMARK. Solution is a vector-function y(t) = (y_1(t), ..., y_n(t)), where n = colCount - 1.
		        The array dataframe has the following structure:
		           dataframe[0, ..., rowCount - 1] contains times, i.e. values {t0, t0 + h, t0 + 2h, ... , t1};
				   dataframe[rowCount, ..., 2 * rowCount - 1] contains values of {y_1(t0), y_1(t0 + h), ..., y_1(t1)};
				   dataframe[2 * rowCount, ..., 3 * rowCount - 1] contains values of {y_2(t0), y_2(t0 + h), ..., y_2(t1)};
				   . . .
				   dataframe[n * rowCount, ..., (n + 1) * rowCount - 1] contains values of {y_n(t0), y_n(t0 + h), ..., y_n(t1)};
	*/
	template<typename DataType, typename ArgType, typename VecType>
	int solveODE(VecType(*f)(ArgType, VecType&), 
		DataType t0, DataType t1, DataType h, DataType * y0, DataType tol,
		DataType * dataFrame, int rowCount, int colCount)
	{
		// dimension of solution
		int dim = colCount - 1;

		// operating variables
		ArgType _t0 = static_cast<ArgType>(t0);
		ArgType _t1 = static_cast<ArgType>(t1);
		ArgType _hInitial = static_cast<ArgType>(h);
		ArgType _tol = static_cast<ArgType>(tol);

		// vector - initial value of y, i.e. y0
		VecType yInitial(dim);
		for (int i = 0; i < dim; i++)
			yInitial(i) = y0[i];

		// structures for {t} and {y(t)}: applying vector provides higher performance than list and deque
		std::vector<ArgType> times;
		std::vector<VecType> solutions;

		/*std::list<ArgType> times;
		std::list<VecType> solutions;*/

		/*std::deque<ArgType> times;
		std::deque<VecType> solutions;*/

		// solve ODE: times and solutions are obtained
		int resultCode = RKCKsolver(f, _t0, _t1, _hInitial, yInitial, _tol, times, solutions);
		if (resultCode != NO_ERRORS)
			return resultCode;

		// interpolation of results
		resultCode = linearInterpolation(_t0, _t1, _hInitial, dataFrame, rowCount, colCount, times, solutions);

		return resultCode;
	} // solve


}; // ode


#endif


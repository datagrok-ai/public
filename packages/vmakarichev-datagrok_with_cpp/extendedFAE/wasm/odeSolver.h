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

#include <vector>

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

namespace ode {

	enum ResultCode {NO_ERRORS = 0, UNKNOWN_PROBLEM, METHOD_FAILS };	

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
		VecType& y, VecType& dydt, ArgType t, ArgType h, VecType& yOut, VecType& yErr ) noexcept
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
		VecType& y, VecType& dydt, ArgType& t, ArgType hTry, VecType& yScale, ArgType& hNext, ArgType tol ) noexcept
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
		ArgStruct & times, VecStruct & solutions) noexcept
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
		bool flag = true;

		// compute numerical solution
		while (flag)
		{
			// compute derivative
			VecType dydt = f(t, y);
			
			// compute scale vector
			VecType yScale = y.cwiseAbs() + (dydt * h).cwiseAbs() + tiny;

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
		ArgStruct & times, VecStruct & solutions) noexcept
	{
		// Copying solution at the point t0
		dataFrame[0] = static_cast<DataType>(times[0]);

		for (int j = 1; j < colCount; j++)
			dataFrame[j * rowCount] = static_cast<DataType>(solutions[0](j - 1));

		int structIndexLimit = times.size() - 1;
		int rowLimit = rowCount - 1;

		// Copying solution at the point t1
		dataFrame[rowLimit] = static_cast<DataType>(times[structIndexLimit]);

		for (int j = 1; j < colCount; j++)
			dataFrame[rowLimit + j * rowCount] = static_cast<DataType>(solutions[structIndexLimit](j - 1));
				
		int m = 1;

		OperatingType t = t0 + h;		

		// the classic linear interpolation method: data is computed for the rest of points		
		for(int i = 1; i < rowLimit; i++)
		{
			while (t > times[m])
				m++;

			OperatingType cLeft = (times[m] - t) / (times[m] - times[m - 1]);
			OperatingType cRight = 1.0 - cLeft;

			dataFrame[i] = static_cast<DataType>(t);

			for (int j = 1; j < colCount; j++)
			{
				dataFrame[i + j * rowCount]
					= static_cast<DataType>(cRight * solutions[m](j - 1)
						+ cLeft * solutions[m - 1](j - 1));
			}

			t += h;
		} // for i

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
		DataType * dataFrame, int rowCount, int colCount) noexcept
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

		// structures for {t} and {y(t)}
		std::vector<ArgType> times;
		std::vector<VecType> solutions;

		times.reserve(rowCount);
		solutions.reserve(rowCount);

		// solve ODE: times and solutions are obtained
		int resultCode = RKCKsolver(f, _t0, _t1, _hInitial, yInitial, _tol, times, solutions);
		if (resultCode != NO_ERRORS)
			return resultCode;

		// interpolation of results
		resultCode = linearInterpolation(_t0, _t1, _hInitial, dataFrame, rowCount, colCount, times, solutions);

		return resultCode;
	} // solveODE

	/*  Solver of the initial ODE problem dy/dt = f(t,y), y(t0) = y0.	    
		Solves the problem on the segment [t0, t1], and linearly interpolated are stored in the dataframe.
		Spatial complexity reduction is applied: additional memory costs are O(1).
		  f - right part of the equation;
		  t0 - initial point;
		  t1 - end point;
		  step - step;
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
	int solveODElight(VecType(*f)(ArgType, VecType&),
		DataType t0, DataType t1, const DataType step, DataType* y0, DataType tol,
		DataType* dataFrame, int rowCount, int colCount)
	{
		// Here, Cash-Karp method is implemented (see [2] for more details).

		// dimension of solution
		int dim = colCount - 1;

		// operating variables
		ArgType _t0 = static_cast<ArgType>(t0);
		ArgType _t1 = static_cast<ArgType>(t1);
		ArgType h = static_cast<ArgType>(step);
		ArgType hDataframe = static_cast<ArgType>(step);
		ArgType tolerance = static_cast<ArgType>(tol);

		// vector - initial value of y, i.e. y0
		VecType yInitial(dim);
		for (int i = 0; i < dim; i++)
			yInitial(i) = y0[i];
				
		// method routine
		VecType tiny = VecType::Ones(dim) * 1e-20;
		ArgType timeDataframe = _t0 + hDataframe;

		ArgType t = _t0;
		ArgType tPrev = _t0;
		ArgType hNext = 0.0;
		VecType y = yInitial;
		VecType yPrev = yInitial;
		VecType dydt(dim);
		VecType yScale(dim);
		bool flag = true;
		int index = 1;

		// 1. solution at the point t0
		dataFrame[0] = t0;

		for (int i = 1; i <= dim; i++)
			dataFrame[i * rowCount] = y0[i - 1];

		// compute numerical solution for the points from the interval (t0, t1)
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
			int resultCode = adaptiveRKCK(f, y, dydt, t, h, yScale, hNext, tolerance);

			// check result of one step
			if (resultCode != NO_ERRORS)
				return resultCode;
						
			//cout << "  " << t << " <--> " << y.transpose() << endl;

			// compute lineraly interpolated results and store them in dataframe
			while (timeDataframe < t)
			{
				ArgType cLeft = (t - timeDataframe) / (t - tPrev);
				ArgType cRight = 1.0 - cLeft;

				dataFrame[index] = static_cast<DataType>(timeDataframe);

				for (int j = 1; j < colCount; j++)
				{
					dataFrame[index + j * rowCount]
						= static_cast<DataType>(cRight * y(j - 1)
							+ cLeft * yPrev(j - 1));
				}

				//cout << timeDataframe << "   " << step << endl;

				timeDataframe += hDataframe;
				index++;
			}
			

			h = hNext;

			tPrev = t;
			yPrev = y;
		} // while	

		// 3. solution at the point t1
		dataFrame[rowCount - 1] = t1;

		for (int i = 1; i <= dim; i++)
			dataFrame[(i + 1) * rowCount - 1] = static_cast<DataType>(y(i - 1));

		return NO_ERRORS;
	} // solveODElight


}; // ode


#endif


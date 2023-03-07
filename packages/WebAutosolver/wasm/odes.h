// odes.h

/* Adaptive step solvers for the initial problems for ODEs:
     - explicit Runge-Kutta Cash-Karp method;
     - implcit modified Rosnbrock triple method.

   The formulas applied are taken from

     [1] Endre Suli and David F. Mayers. An Introduction to Numerical Analysis, 2003.

     [2] Steve Chapra and Raymond P. Canale. Numerical Methods for Engineers, 2021.

     [3] https://doi.org/10.1137/S1064827594276424

     [4] https://doi.org/10.1016/S0898-1221(00)00175-9
*/

#ifndef ODES_H
#define ODES_H

#include "../../../../../Eigen/Eigen/Dense"
using namespace Eigen;

namespace odes {

	enum ResultCode { NO_ERRORS = 0, UNKNOWN_PROBLEM, METHOD_FAILS };

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
	int RKCK(VecType(*f)(ArgType, VecType&),
		VecType& y, VecType& dydt, ArgType t, ArgType h, VecType& yOut, VecType& yErr)
	{
		// implementation of R.-K. Cash-Karp method (see [2] for more details)

		//TODO: move method constants to separate file
		//cout << "\ny = " << y.transpose() << endl;
		VecType ytemp = y + 0.2 * h * dydt;
		//cout << "\nytemp = " << ytemp.transpose() << endl;
		VecType k2 = f(t + 0.2 * h, ytemp);
		//cout << "\nk2 = " << k2.transpose() << endl;
		ytemp = y + h * (0.075 * dydt + 0.225 * k2);
		//cout << "\nytemp = " << ytemp.transpose() << endl;
		VecType k3 = f(t + 0.3 * h, ytemp);
		//cout << "\nk3 = " << k3.transpose() << endl;
		ytemp = y + h * (0.3 * dydt - 0.9 * k2 + 1.2 * k3);
		//cout << "\nytemp = " << ytemp.transpose() << endl;
		VecType k4 = f(t + 0.6 * h, ytemp);
		//cout << "\nk4 = " << k4.transpose() << endl;
		ytemp = y + h * (-0.2037037037037037 * dydt + 2.5 * k2 - 2.592592592592593 * k3 + 1.296296296296296 * k4);
		//cout << "\nytemp = " << ytemp.transpose() << endl;
		VecType k5 = f(t + 0.6 * h, ytemp);
		ytemp = y + h * (0.0294958043981481 * dydt + 0.341796875 * k2 + 0.0415943287037037 * k3
			+ 0.4003454137731481 * k4 + 0.061767578125 * k5);
		VecType k6 = f(t + 0.875 * h, ytemp);
		yOut = y + h * (0.0978835978835979 * dydt + 0.4025764895330113 * k3
			+ 0.21043771043771045 * k4 + 0.2891022021456804 * k6);
		yErr = h * (-0.004293774801587311 * dydt + 0.018668586093857853 * k3 - 0.034155026830808066 * k4
			- 0.019321986607142856 * k5 + 0.03910220214568039 * k6);

		//cin.get();

		return NO_ERRORS;
	} // RKCK


	/*  One stage of modifed Rosenbrock triple method for solving the problem dy/dt = f(t,y), y(t0) = y0.
		Computes yOut and yErr, which are respectively numerical solution and errors at the point t + h.
		  f - right part of the ODE solved;
		  y - solution at the point t;
		  dydt - derivative of the solution at the point t;
		  t - current value of the independent variable;
		  h - step;
		  yOut - numerical solution at the point t + h;
		  yErr - errors at the point t + h.
	*/
	template<typename ArgType, typename VecType, typename MatType>
	int MRT(VecType(*f)(ArgType, VecType&),
		VecType(*T)(ArgType, VecType&, ArgType),
		MatType(*J)(ArgType, VecType&, ArgType),
		VecType& y, VecType& dydt, ArgType t, ArgType h, VecType& yOut, VecType& yErr)
	{
		const ArgType EPS = 1.0e-10;// 0.000001;

		auto dimCount = y.size();

		// Quantities used in Rosenbrock method (see papers [3, 4] for more details)
		const ArgType d = static_cast<ArgType>(1.0 - sqrt(2.0) / 2.0);
		const ArgType e32 = static_cast<ArgType>(6.0 + sqrt(2.0));
		MatType I = MatType::Identity(dimCount, dimCount);
		MatType W(dimCount, dimCount);
		MatType invW(dimCount, dimCount);
		VecType f0(dimCount);
		VecType k1(dimCount);
		VecType f1(dimCount);
		VecType k2(dimCount);
		VecType f2(dimCount);
		VecType k3(dimCount);
		VecType yDer(dimCount);
		MatType hdT = h * d * T(t, y, EPS);

		// The main computations
		f0 = f(t, y);
		W = I - h * d * J(t, y, EPS);
		invW = W.inverse();
		k1 = invW * (f0 + hdT); //h * d * T(t, y, EPS));
		yDer = y + 0.5 * h * k1;
		f1 = f(t + 0.5 * h, yDer);
		k2 = invW * (f1 - k1) + k1;
		yOut = y + k2 * h;
		f2 = f(t + h, yOut);
		k3 = invW * (f2 - e32 * (k2 - f1) - 2.0 * (k1 - f0) + hdT); // h * d * T(t, y, EPS));
		yErr = (k1 - 2.0 * k2 + k3) * h / 6;

		return NO_ERRORS;
	} // MRT

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
	template<typename ArgType, typename VecType, typename MatType>
	int adaptiveMRT(VecType(*f)(ArgType, VecType&),
		VecType(*T)(ArgType, VecType&, ArgType),
		MatType(*J)(ArgType, VecType&, ArgType),
		VecType& y, VecType& dydt, ArgType& t, ArgType hTry, VecType& yScale, ArgType& hNext, ArgType tol)
	{
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
			// one stage of the modified Rosenbrok triple approach
			int resultCode = MRT(f, T, J, y, dydt, t, h, yTemp, yErr);

			// check result code
			if (resultCode != NO_ERRORS)
				return resultCode;

			// estimating error
			ArgType errmax = 0.0;
			for (unsigned i = 0; i < dim; i++)
				errmax = fmax(errmax, fabs(yErr(i) / yScale(i)));
			errmax /= tol;

			//cout << "errmax: " << errmax << endl;

			/*t += h;
			hNext = h;
			break;*/

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
	} // adaptiveMRT

	/*  Solver of the initial ODE problem dy/dt = f(t,y), y(t0) = y0.
		Solves the problem on the segment [t0, t1], and linearly interpolated are stored in the dataframe.
		Implicit adaptive step modified Rosenbrok triple method is applied.
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
	template<typename DataType, typename ArgType, typename VecType, typename MatType>
	int solveODE(VecType(*f)(ArgType, VecType&),
		VecType(*T)(ArgType, VecType&, ArgType),
		MatType(*J)(ArgType, VecType&, ArgType),
		DataType t0, DataType t1, const DataType step, DataType* y0, DataType tol,
		DataType* dataFrame, int rowCount, int colCount)
	{
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

		unsigned counter = 0;

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

			// call of adaptive step modified Rosenbrok triple method
			int resultCode = adaptiveMRT(f, T, J, y, dydt, t, h, yScale, hNext, tolerance);

			counter++;

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
			} // while

			h = hNext;

			tPrev = t;
			yPrev = y;
		} // while

		// 3. solution at the point t1
		dataFrame[rowCount - 1] = t1;

		for (int i = 1; i <= dim; i++)
			dataFrame[(i + 1) * rowCount - 1] = static_cast<DataType>(y(i - 1));

		//cout << "\nNumber of points found by adaptive method: " << counter << endl;

		return NO_ERRORS;
	} // solveODE

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
		VecType& y, VecType& dydt, ArgType& t, ArgType hTry, VecType& yScale, ArgType& hNext, ArgType tol)
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
		auto dim = y.size();
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


	/*  Solver of the initial ODE problem dy/dt = f(t,y), y(t0) = y0.
		Solves the problem on the segment [t0, t1], and linearly interpolated are stored in the dataframe.
		Explicit adaptive step Runge-Kutta Cash-Karp method is applied.
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
	int solveODE(VecType(*f)(ArgType, VecType&),
		DataType t0, DataType t1, const DataType step, DataType* y0, DataType tol,
		DataType* dataFrame, int rowCount, int colCount)
	{
		// Here, Cash-Karp method is implemented (see [2] for more details).

		// dimension of solution
		int dim = colCount - 1;

		// operating variables
		ArgType _t0 = static_cast<ArgType>(t0);
		ArgType _t1 = static_cast<ArgType>(t1);
		ArgType h = 0.000001;//static_cast<ArgType>(step);
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

		unsigned counter = 0;

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

			counter++;

			// check result of one step
			if (resultCode != NO_ERRORS)
				return resultCode;

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

		//cout << "\nNumber of points found by adaptive method: " << counter << endl;

		return NO_ERRORS;
	} // solveODE


}; // ode


#endif

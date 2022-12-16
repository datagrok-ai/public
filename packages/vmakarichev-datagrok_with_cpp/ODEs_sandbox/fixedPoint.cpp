#include<vector>
#include<iostream>
using namespace std;

/* Computation of the coefficients k_{i,j}, required for non-explicit Runge-Kutta method.
     func - right part of the system of ODEs solved, i.e. dy/dt = func(t,y);
	 t - value of t[n], previous time;
	 y - values of y(t[n]) <- these values are computed at the previous step of R.-K. method;	
	 k - coefficients { k_{i,j} } that are computed;
	 a, c - parameters of R.-K. method;
	 h - step (h = t[n + 1] - t[n]);
	 precision - precision of the applied iteration approach.
*/
void getRoots(vector<double>(*func) (double, vector<double> &),
	double t,
	vector<double> & y,
	vector<vector<double>> & k,
	vector<vector<double>> & a,
	vector<double> & c,
	double h,
	double precision)
{
	int N = y.size(); // number of equations
	int s = c.size(); // number coefficients k-s for each y_i;

	vector<vector<double>> kPrev(N);

	// fill kPrev with zeros
	for (int i = 0; i < N; i++)
	{
		kPrev[i].resize(s);
		for (int j = 0; j < s; j++)
			kPrev[i][j] = 0.0;
	}

	// input to the func
	double tArg = 0; // value of the independant variable
	vector<double> args(N); // dependent variables
	vector<double> argsDer(N);

	double hDer = 0.0001;

	int iter = 1;

	// the main iterations
	while (true)
	{
		// compute each column of the matrix k
		for (int i = 0; i < s; i++)
		{
			tArg = t + c[i] * h;

			// compute args
			for (int r = 0; r < N; r++)
			{
				double sum = 0.0;
				
				for (int j = 0; j < s; j++)
					sum += a[i][j] * kPrev[r][j];

				args[r] = y[r] + sum * h;
			}

			auto res = func(tArg, args);

			for (int r = 0; r < N; r++)
			{
				double sum = 0.0;

				for (int j = 0; j < s; j++)
					sum += a[i][j] * (kPrev[r][j] + hDer);

				argsDer[r] = y[r] + sum * h;
			}
			auto resDer = func(tArg, argsDer);

			// copy res to the column K
			for (int m = 0; m < N; m++) {
				//k[m][i] = res[m] ;
				double function = k[m][i] - res[m];
				double functionAdd = k[m][i] + hDer - resDer[m];
				k[m][i] = k[m][i] - function * hDer / (functionAdd - function);

				/*double dif = k[m][i] - res[m];

				if (k[m][i] != 0.0)
				{
					dif = fabs(dif / k[m][i]);

					if(dif < 1)
						k[m][i] += 0.1 * res[m];
					else
						k[m][i] += 0.001 * res[m];
				}
				else
					k[m][i] += 0.001 * res[m];		*/		
			}
		}		

		// compute maximum absolute deviation between current and previous k-s
		double mad = fabs(k[0][0] - kPrev[0][0]);

		for (int i = 0; i < N; i++)
			for (int j = 0; j < s; j++)
				mad = fmax(mad, k[i][j] - kPrev[i][j]);

		cout << " iteration no. " << iter << "  mad = " << mad << endl;
		iter++;

		// check precision
		if (mad < precision)
			break;
		
		// copy values of k to kPrev
		for (int i = 0; i < N; i++)
			for (int j = 0; j < s; j++)
				kPrev[i][j] = k[i][j];
	} // while
} // getRoots
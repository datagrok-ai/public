// wrapper.cpp
// Declaration of ODEs solver wrapper

#include<vector>
using std::vector;

#include "solver.h"

int solverWrapper(float * tVals, int tValsLength,
	float * params, int paramsLength,
	float * solution, int solutionRowCount, int solutionColumnCount)
{
	vector<double> times(tValsLength);

	for (int i = 0; i < tValsLength; i++)
		times[i] = tVals[i];

	vector<double> _parameters(paramsLength);
	
	for (int i = 0; i < paramsLength; i++)
		_parameters[i] = params[i];

	vector<vector<double>> result;

	Solve(&result, &times, &_parameters, NULL, NULL, NULL);

	for (int i = 0; i < solutionColumnCount; i++)
		for (int j = 0; j < solutionRowCount; j++)
			solution[i * solutionRowCount + j] = static_cast<float>(result[i][j]);
	
	return 0;
}

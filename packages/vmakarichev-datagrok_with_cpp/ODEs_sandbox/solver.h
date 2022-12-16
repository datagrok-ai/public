// solver.h
// Declaration of the function for solving system of ODEs.

#ifndef SOLVER_H
#define SOLVER_H

#include<vector>

void Solve(
	std::vector<std::vector<double>> * result,
	std::vector<double> * times,
	std::vector<double> * _parameters,
	std::vector<double> * _doseTime,
	std::vector<int> * _doseADM,
	std::vector<double> * _doseAMT);

#endif // !SOLVER_H


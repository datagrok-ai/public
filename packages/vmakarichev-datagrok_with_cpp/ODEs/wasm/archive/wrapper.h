// wrapper.h
// Declaration of ODEs solver wrapper

#ifndef WRAPPER_H
#define WRAPPER_H

/* Wrapper for ODEs solver.     
	 tVals - corresponds to solver times;
	 tValsLength - size of tVals;
	 params - corresponds to solver _parameters;
	 paramsLength - size of params;
	 solution - array that contains solver result (result data is concatenated), data for DATAGROK columnList;
	 solutionRowCount - number of rows;
	 solutionColumnCount - number of columns;
*/
int solverWrapper(float * tVals, int tValsLength, 
	float * params, int paramsLength,
	float * solution, int solutionRowCount, int solutionColumnCount);

#endif

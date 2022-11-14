#include<iostream>
using namespace std;

//#include "../../../Eigen/Eigen/Dense"
//using namespace Eigen;

#include "tests.h"

int main()
{
	//differentTypeValuesManipulation();

	//pcaUsingCorrelationMatrixTest();

	//comparePerformanceOfPCAwithCorMatr();

	//investigatePCAwithCorMatrVoidDataRepresentation(100);

	investigatePCAwithCorMatrFloatDataRepresentation(10);

	cin.get();

	return 0;
}
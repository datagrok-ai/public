#include<iostream>
using namespace std;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "tests.h"

int main()
{
	//testCreateDataset();

	//testTrainModelSimpleLinear();

	//testGeneratorLinearSeparable();

	//testGeneratorLinearNonSeparable();

	testTrainModelComplexLinear();

	return 0;
}
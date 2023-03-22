#include<iostream>
using namespace std;

#include "../../../Eigen/Eigen/Dense"
using namespace Eigen;

#include "tests.h"

int main()
{
	//testCreateDataset();

	//testCreateNormalizedDataset();

	//testTrainModelSimpleLinear();

	//testGeneratorLinearSeparable();

	//testGeneratorLinearNonSeparable();

	//testTrainModelNormalizedDataLinear();

	testTrainModelNormalizedDataLinearHighDim(false);

	//testTrainModelNormalizedDataLinearHighDimDouble();

	return 0;
}
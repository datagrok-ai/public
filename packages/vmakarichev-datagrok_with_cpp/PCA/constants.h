#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace pca
{
	// Constants for performance investigation
	const int SEED = 1234;
	const int WIDTH = 100000;
	const int HEIGHT_OF_FLOATS = 50;
	const int HEIGHT_OF_INTS = 50;
	const int HEIGHT = HEIGHT_OF_FLOATS + HEIGHT_OF_INTS;
	const int NUM_OF_PRINCIPAL_COMPONENTS = HEIGHT;
};

#endif // ! CONSTANTS_H


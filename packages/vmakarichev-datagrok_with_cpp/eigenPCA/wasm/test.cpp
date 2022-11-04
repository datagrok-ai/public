#include<iostream>
using namespace std;

#include "pcaExport.cpp"

int main()
{
    cout << "Hi!";

    const int HEIGHT = 5;
	const int WIDTH = 4;
	const int SIZE = HEIGHT * WIDTH;
	const int NUM_OF_PRINCIPAL_COMPONENTS = HEIGHT;
	const int PRINCIPAL_COMPONENTS_SIZE = NUM_OF_PRINCIPAL_COMPONENTS * WIDTH;

	/*float data[SIZE] = { 11, 12, 13, 14, 15,
						16, 17, 18, 19, 20,
						21, 22, 23, 24, 25,
						26, 27, 28, 29, 30 };*/

	float data[SIZE] = { 141, 212, -313, 314, -215,
		                 166, -717, -118, 19, 120,
		                 -321, 202, 123, -24, 325,
		                 276, -227, 218, -129, 130 };

	cout << "data:\n";
	for (int i = 0; i < SIZE; i++)
		cout << "  " << data[i];

	float princapalComponents[PRINCIPAL_COMPONENTS_SIZE] = { 0 };

	cout << "\n\nprincipal components (before):\n";
	for (int i = 0; i < PRINCIPAL_COMPONENTS_SIZE; i++)
		cout << "  " << princapalComponents[i];

	float approximation[SIZE] = { 0 };

	cout << "\n\napproximation (before):\n";
	for (int i = 0; i < SIZE; i++)
		cout << "  " << approximation[i];

	cout << "\n\nPCA ...\n";

	principalComponentAnalysis(data, HEIGHT, WIDTH, NUM_OF_PRINCIPAL_COMPONENTS, princapalComponents, approximation);

	cout << "\n\nprincipal components (after):\n";
	for (int i = 0; i < PRINCIPAL_COMPONENTS_SIZE; i++)
		cout << "  " << princapalComponents[i];
	
	cout << "\n\napproximation (after):\n";
	for (int i = 0; i < SIZE; i++)
		cout << "  " << approximation[i];

	//cout << "\n\nmaximum absolute deviation: " << mad(data, approximation, SIZE);

    return 0;
}
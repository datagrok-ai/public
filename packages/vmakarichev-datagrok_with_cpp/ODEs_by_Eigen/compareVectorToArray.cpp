#include<iostream>
#include<vector>
#include<ctime>
#include<cstdlib>
using namespace std;

// vector vs pointer array
void compareVectorToArray()
{
	unsigned N = 2;

	double* arr = new double[N];

	vector<double> vec(N);

	for (unsigned i = 0; i < N; i++)
		arr[i] = vec[i] = rand() % 10;

	auto start = time(0);

	double sum = 0;

	for (unsigned i = 0; i < 10000000; i++)
		for (unsigned j = 0; j < 10000000; j++)
			for (unsigned k = 0; k < 10000000; k++)
				for (unsigned m = 0; m < 10000000; m++)
					//sum = arr[0] + arr[1];
					sum = vec[0] + vec[1];
					//sum = vec.at(0) + vec.at(1);

	auto finish = time(0);

	cout << "\ntime is " << finish - start << " sec.\n";

	delete[] arr;
}
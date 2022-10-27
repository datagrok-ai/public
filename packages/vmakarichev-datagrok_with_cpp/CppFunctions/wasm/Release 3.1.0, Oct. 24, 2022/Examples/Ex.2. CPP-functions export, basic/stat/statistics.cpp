#include "statistics.h"

float sta::minOfArray(const float * arr, const int arrLength)
{
	float min = arr[0];

	for (int i = 0; i < arrLength; ++i)
		if (min > arr[i])
			min = arr[i];

	return min;
}

float sta::maxOfArray(const float * arr, const int arrLength)
{
	float max = arr[0];

	for (int i = 0; i < arrLength; ++i)
		if (max < arr[i])
			max = arr[i];

	return max;
}
#ifndef INDICATORS_H
#define INDICATORS_H

namespace ind {

	template <class T>
	T minOfArray(const T * arr, const int length)
	{
		T min = arr[0];

		for (int i = 1; i < length; ++i)
			if (min > arr[i])
				min = arr[i];

		return min;
	}

	template <class T>
	T maxOfArray(const T * arr, const int length)
	{
		T max = arr[0];

		for (int i = 1; i < length; ++i)
			if (max < arr[i])
				max = arr[i];

		return max;
	}
};

#endif 

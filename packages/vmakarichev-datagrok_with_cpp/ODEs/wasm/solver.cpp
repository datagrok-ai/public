// solver.h
// Implementation of the function for solving system of ODEs.

#include <cmath>
#include <string>
#include <vector>

#include "solver.h"

void Solve(
	std::vector<std::vector<double>> * result,
	std::vector<double> * times,
	std::vector<double> * _parameters,
	std::vector<double> * _doseTime,
	std::vector<int> * _doseADM,
	std::vector<double> * _doseAMT)
{
	int points = times->size();
	std::vector<double>::iterator _itB_times = times->begin();
	std::vector<double>::iterator _itE_times = times->end();
	double _time = *_itB_times;
	double step = 0;

	//parameters
	double VLinitial = (*_parameters)[0];
	//constants
	double k1red = 0.000934 * 60;
	double k1ox = 0.00018 * 60;
	double MEAthiol = 34;
	double pH = 7.4;
	double pKa2MEA = 8.19;
	//allocation for selected variables
	std::vector<double> * FFred = new std::vector<double>();
	std::vector<double> * KKred = new std::vector<double>();
	std::vector<double> * FFox = new std::vector<double>();
	std::vector<double> * KKox = new std::vector<double>();
	std::vector<double> * VL = new std::vector<double>();
	//initial expressions and starting values
	double MEAthiolate_t = MEAthiol*std::pow(10, (pH - pKa2MEA));
	double FFred_t = 0;
	double KKred_t = 0;
	double FFox_t = 50;
	double KKox_t = 50;
	double VL_t = 6.6;
	//selected variables to memory
	FFred->push_back(FFred_t);
	KKred->push_back(KKred_t);
	FFox->push_back(FFox_t);
	KKox->push_back(KKox_t);
	VL->push_back(VL_t);
	++_itB_times;

	//dosages at 0 time point
	int _event = 0;
	//ODE Solver: explicit Sofroniou-Spaletta 4 order method
	double FFred_k1 = 0;
	double KKred_k1 = 0;
	double FFox_k1 = 0;
	double KKox_k1 = 0;
	double VL_k1 = 0;
	double FFred_k2 = 0;
	double KKred_k2 = 0;
	double FFox_k2 = 0;
	double KKox_k2 = 0;
	double VL_k2 = 0;
	double FFred_k3 = 0;
	double KKred_k3 = 0;
	double FFox_k3 = 0;
	double KKox_k3 = 0;
	double VL_k3 = 0;
	double FFred_add = 0;
	double KKred_add = 0;
	double FFox_add = 0;
	double KKox_add = 0;
	double VL_add = 0;
	double MEAthiolate_k1 = 0;
	double MEAthiolate_k2 = 0;
	double MEAthiolate_k3 = 0;
	double MEAthiolate_add = 0;

	for (_itB_times; _itB_times != _itE_times; ++_itB_times)
	{
		step = *_itB_times - _time;
		_time = *_itB_times;

		FFred_k1 = 0;
		KKred_k1 = 0;
		FFox_k1 = -k1red*FFox_t*(VLinitial / VL_t)*std::pow(MEAthiolate_t, 2)*std::pow((VLinitial / VL_t), 2) + k1ox*FFred_t*(VLinitial / VL_t);
		KKox_k1 = -k1red*KKox_t*(VLinitial / VL_t)*std::pow(MEAthiolate_t, 2)*std::pow((VLinitial / VL_t), 2) + k1ox*KKred_t*(VLinitial / VL_t);
		VL_k1 = 0;
		FFred_add = FFred_t + 0.5*step*FFred_k1;
		KKred_add = KKred_t + 0.5*step*KKred_k1;
		FFox_add = FFox_t + 0.5*step*FFox_k1;
		KKox_add = KKox_t + 0.5*step*KKox_k1;
		VL_add = VL_t + 0.5*step*VL_k1;
		MEAthiolate_add = MEAthiol*std::pow(10, (pH - pKa2MEA));

		FFred_k2 = 0;
		KKred_k2 = 0;
		FFox_k2 = -k1red*FFox_add*(VLinitial / VL_add)*std::pow(MEAthiolate_add, 2)*std::pow((VLinitial / VL_add), 2) + k1ox*FFred_add*(VLinitial / VL_add);
		KKox_k2 = -k1red*KKox_add*(VLinitial / VL_add)*std::pow(MEAthiolate_add, 2)*std::pow((VLinitial / VL_add), 2) + k1ox*KKred_add*(VLinitial / VL_add);
		VL_k2 = 0;
		FFred_add = FFred_t + 0.5*step*FFred_k2;
		KKred_add = KKred_t + 0.5*step*KKred_k2;
		FFox_add = FFox_t + 0.5*step*FFox_k2;
		KKox_add = KKox_t + 0.5*step*KKox_k2;
		VL_add = VL_t + 0.5*step*VL_k2;
		MEAthiolate_add = MEAthiol*std::pow(10, (pH - pKa2MEA));

		FFred_k3 = 0;
		KKred_k3 = 0;
		FFox_k3 = -k1red*FFox_add*(VLinitial / VL_add)*std::pow(MEAthiolate_add, 2)*std::pow((VLinitial / VL_add), 2) + k1ox*FFred_add*(VLinitial / VL_add);
		KKox_k3 = -k1red*KKox_add*(VLinitial / VL_add)*std::pow(MEAthiolate_add, 2)*std::pow((VLinitial / VL_add), 2) + k1ox*KKred_add*(VLinitial / VL_add);
		VL_k3 = 0;
		FFred_add = FFred_t + step*FFred_k3;
		KKred_add = KKred_t + step*KKred_k3;
		FFox_add = FFox_t + step*FFox_k3;
		KKox_add = KKox_t + step*KKox_k3;
		VL_add = VL_t + step*VL_k3;
		MEAthiolate_add = MEAthiol*std::pow(10, (pH - pKa2MEA));

		FFred_t += step*(FFred_k1 + 2 * FFred_k2 + 2 * FFred_k3 + (0)) / 6;
		KKred_t += step*(KKred_k1 + 2 * KKred_k2 + 2 * KKred_k3 + (0)) / 6;
		FFox_t += step*(FFox_k1 + 2 * FFox_k2 + 2 * FFox_k3 + (-k1red*FFox_add*(VLinitial / VL_add)*std::pow(MEAthiolate_add, 2)*std::pow((VLinitial / VL_add), 2) + k1ox*FFred_add*(VLinitial / VL_add))) / 6;
		KKox_t += step*(KKox_k1 + 2 * KKox_k2 + 2 * KKox_k3 + (-k1red*KKox_add*(VLinitial / VL_add)*std::pow(MEAthiolate_add, 2)*std::pow((VLinitial / VL_add), 2) + k1ox*KKred_add*(VLinitial / VL_add))) / 6;
		VL_t += step*(VL_k1 + 2 * VL_k2 + 2 * VL_k3 + (0)) / 6;

		MEAthiolate_t = MEAthiol*std::pow(10, (pH - pKa2MEA));

		FFred->push_back(FFred_t);
		KKred->push_back(KKred_t);
		FFox->push_back(FFox_t);
		KKox->push_back(KKox_t);
		VL->push_back(VL_t);
	}

	//selected variables to memory
	result->push_back(*FFred);
	result->push_back(*KKred);
	result->push_back(*FFox);
	result->push_back(*KKox);
	result->push_back(*VL);
}


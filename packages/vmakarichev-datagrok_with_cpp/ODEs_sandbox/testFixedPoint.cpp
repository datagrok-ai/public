#include<iostream>
using namespace std;

#include "test.h"
#include "fixedPoint.h"

vector<double> func(double t, vector<double> & y)
{
	vector<double> res(y.size());

	res[0] = y[0] + y[1];
	res[1] = y[0] - y[1];

	return res;
}

// test computation of k-s using fixedPoint
void testFixedPoint()
{
	int N = 2;
	int s = 4;

	double t = 1.3;

	vector<double> y(N);
	y[0] = 2.1;
	y[1] = 1.2;

	vector<vector<double>> k(N);
	for (int i = 0; i < N; i++)
		k[i].resize(s);

	cout << "\nk:\n";
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < s; j++)
			cout << "  " << k[i][j];
		cout << endl;
	}

	vector<vector<double>> a(s);
	for (int i = 0; i < s; i++)
		a[i].resize(s);

	a[0][0] = 0.1; a[0][1] = 0.2; a[0][2] = 0.3; a[0][3] = 0.2;
	a[1][0] = 0.2; a[1][1] = 0.3; a[1][2] = 0.2; a[1][3] = 0.1;
	a[2][0] = 0.1; a[2][1] = 0.1; a[2][2] = 0.2; a[2][3] = 0.1;
	a[3][0] = 0.1; a[3][1] = 0.2; a[3][2] = 0.1; a[3][3] = 0.1;

	cout << "\na:\n";
	for (int i = 0; i < s; i++)
	{
		for (int j = 0; j < s; j++)
			cout << "  " << a[i][j];
		cout << endl;
	}

	vector<double> c(s);
	for (int i = 0; i < s; i++)
	{
		c[i] = a[i][0];

		for (int j = 1; j <= i; j++)
			c[i] += a[i][j];
	}

	cout << "\nc:\n";
	for (int i = 0; i < s; i++)
		cout << "  " << c[i];
	cout << endl;

	double h = 0.1;

	double precision = 0.000001;

	getRoots(func, t, y, k, a, c, h, precision);

	cout << "\nFinal result k:\n";
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < s; j++)
			cout << "  " << k[i][j];
		cout << endl;
	}
}


/*std::vector<double> DerivFunc(double _time, std::vector<double> & _y) {
	std::vector<double> _res(_y.size());
	//constants
	double qin = 1;
	double percentO2saturation = 100;
	double yO2in = 0.209;
	double pKa2MEA = 8.19;
	double H = 1.072;
	double T = 300;
	double R = 0.082;
	double P = 1;
	double VLinitial = 6.6;
	double VL = 6.6;
	double Vtotalvessel = 10;
	double AgitatorSpeed = 200;
	double AgitatorDiameter = 2.36;
	double AgitatorPowerNumber = 2.1;
	double pH = 7.4;
	double k1red = 0.000934 * 60;
	double k1ox = 0.00018 * 60;
	double k2Fd = 0.0225 * 60;
	double k2Fa = 184 * 60;
	double k2Kd = 0.000673 * 60;
	double k2Ka = 200 * 60;
	double k3FKa = 302 * 60;
	double k3FKd = 0.000198 * 60;
	double k4ox = 0.00018 * 60;
	double k4red = 0.000934 * 60;
	double kthiolox = 0.04;
	double krcyst = 0;
	double FFox_t = _y[0];
	double KKox_t = _y[1];
	double FFred_t = _y[2];
	double KKred_t = _y[3];
	double Ffree_t = _y[4];
	double Kfree_t = _y[5];
	double MEAthiol_t = _y[6];
	//expressions
	double MEAthiolate_t = MEAthiol_t*std::pow(10, (pH - pKa2MEA));
	double Fin_t = _time<120 ? (0) : (0.025);
	double Fpermeate_t = _time<120 ? (0.025) : (Fin_t);
	double FFox_to_FFred_t = k1red*FFox_t*(VLinitial / VL)*std::pow(MEAthiolate_t, 2)*std::pow((VLinitial / VL), 2);
	double FFred_to_FFox_t = k1ox*FFred_t*(VLinitial / VL);
	double KKox_to_KKred_t = k1red*KKox_t*(VLinitial / VL)*std::pow(MEAthiolate_t, 2)*std::pow((VLinitial / VL), 2);
	double KKred_to_KKox_t = k1ox*KKred_t*(VLinitial / VL);
	double FFred_to_Ffree_t = k2Fd*FFred_t*(VLinitial / VL);
	double Ffree_to_FFred_t = k2Fa*std::pow(Ffree_t, 2)*std::pow((VLinitial / VL), 2)*std::pow(MEAthiolate_t, 2)*std::pow((VLinitial / VL), 2);
	double KKred_to_Kfree_t = k2Kd*KKred_t*(VLinitial / VL);
	double Kfree_to_KKred_t = k2Ka*std::pow(Kfree_t, 2)*std::pow((VLinitial / VL), 2)*std::pow(MEAthiolate_t, 2)*std::pow((VLinitial / VL), 2);
	_res[0] = -FFox_to_FFred_t + FFred_to_FFox_t;
	_res[1] = -KKox_to_KKred_t + KKred_to_KKox_t;
	_res[2] = FFox_to_FFred_t - FFred_to_FFox_t - FFred_to_Ffree_t + Ffree_to_FFred_t;
	_res[3] = KKox_to_KKred_t - KKred_to_KKox_t - KKred_to_Kfree_t + Kfree_to_KKred_t;
	_res[4] = FFred_to_Ffree_t - Ffree_to_FFred_t;
	_res[5] = KKred_to_Kfree_t - Kfree_to_KKred_t;
	_res[6] = -2 * FFox_to_FFred_t + 2 * FFred_to_FFox_t;
	return _res;
}*/

std::vector<double> DerivFunc(double _time, std::vector<double> & _y) {
	std::vector<double> _res(_y.size());
	//constants
	double qin = 1;
	double percentO2saturation = 100;
	double yO2in = 0.209;
	double pKa2MEA = 8.19;
	double H = 1.072;
	double T = 300;
	double R = 0.082;
	double P = 1;
	double VLinitial = 6.6;
	double VL = 6.6;
	double Vtotalvessel = 10;
	double AgitatorSpeed = 200;
	double AgitatorDiameter = 2.36;
	double AgitatorPowerNumber = 2.1;
	double pH = 7.4;
	double k1red = 0.000934 * 60;
	double k1ox = 0.00018 * 60;
	double k2Fd = 0.0225 * 60;
	double k2Fa = 1840000 * 60;
	double k2Kd = 0.000673 * 60;
	double k2Ka = 2000000 * 60;
	double k3FKa = 3020000 * 60;
	double k3FKd = 0.000198 * 60;
	double k4ox = 0.00018 * 60;
	double k4red = 0.000934 * 60;
	double kthiolox = 0.04;
	double krcyst = 0;
	double FFox_t = _y[0];
	double KKox_t = _y[1];
	double FFred_t = _y[2];
	double KKred_t = _y[3];
	double Ffree_t = _y[4];
	double Kfree_t = _y[5];
	double MEAthiol_t = _y[6];
	//expressions
	double MEAthiolate_t = MEAthiol_t*std::pow(10, (pH - pKa2MEA));
	double Fin_t = _time<120 ? (0) : (0.025);
	double Fpermeate_t = _time<120 ? (0.025) : (Fin_t);
	double FFox_to_FFred_t = k1red*FFox_t*(VLinitial / VL)*std::pow(MEAthiolate_t, 2)*std::pow((VLinitial / VL), 2);
	double FFred_to_FFox_t = k1ox*FFred_t*(VLinitial / VL);
	double KKox_to_KKred_t = k1red*KKox_t*(VLinitial / VL)*std::pow(MEAthiolate_t, 2)*std::pow((VLinitial / VL), 2);
	double KKred_to_KKox_t = k1ox*KKred_t*(VLinitial / VL);
	double FFred_to_Ffree_t = k2Fd*FFred_t*(VLinitial / VL);
	double Ffree_to_FFred_t = k2Fa*std::pow(Ffree_t, 2)*std::pow((VLinitial / VL), 2)*std::pow(MEAthiolate_t, 2)*std::pow((VLinitial / VL), 2);
	double KKred_to_Kfree_t = k2Kd*KKred_t*(VLinitial / VL);
	double Kfree_to_KKred_t = k2Ka*std::pow(Kfree_t, 2)*std::pow((VLinitial / VL), 2)*std::pow(MEAthiolate_t, 2)*std::pow((VLinitial / VL), 2);
	_res[0] = -FFox_to_FFred_t + FFred_to_FFox_t;
	_res[1] = -KKox_to_KKred_t + KKred_to_KKox_t;
	_res[2] = FFox_to_FFred_t - FFred_to_FFox_t - FFred_to_Ffree_t + Ffree_to_FFred_t;
	_res[3] = KKox_to_KKred_t - KKred_to_KKox_t - KKred_to_Kfree_t + Kfree_to_KKred_t;
	_res[4] = FFred_to_Ffree_t - Ffree_to_FFred_t;
	_res[5] = KKred_to_Kfree_t - Kfree_to_KKred_t;
	_res[6] = -2 * FFox_to_FFred_t + 2 * FFred_to_FFox_t - 2 * KKox_to_KKred_t + 2 * KKred_to_KKox_t + 2 * FFred_to_Ffree_t + 2 * KKred_to_Kfree_t - 2 * Ffree_to_FFred_t - 2 * Kfree_to_KKred_t;
	return _res;
}

// test computation of k-s using fixedPoint, another example
void testFixedPointAnother()
{
	int N = 7;
	int s = 1;

	double t = 2.331;

	vector<double> y(N);
	y[0] = 2.1;
	y[1] = 1.2;
	y[2] = 1.2;
	y[3] = 2.1;
	y[4] = 1.2;
	y[5] = 1.2;
	y[6] = 1.2;

	vector<vector<double>> k(N);
	for (int i = 0; i < N; i++)
		k[i].resize(s);

	cout << "\nk:\n";
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < s; j++)
			cout << "  " << k[i][j];
		cout << endl;
	}

	vector<vector<double>> a(s);
	for (int i = 0; i < s; i++)
		a[i].resize(s);

	a[0][0] = 1;

	cout << "\na:\n";
	for (int i = 0; i < s; i++)
	{
		for (int j = 0; j < s; j++)
			cout << "  " << a[i][j];
		cout << endl;
	}

	vector<double> c(s);
	for (int i = 0; i < s; i++)
	{
		c[i] = a[i][0];

		for (int j = 1; j <= i; j++)
			c[i] += a[i][j];
	}

	cout << "\nc:\n";
	for (int i = 0; i < s; i++)
		cout << "  " << c[i];
	cout << endl;

	double h = 0.00001;

	double precision = 0.01;

	getRoots(DerivFunc, t, y, k, a, c, h, precision);

	cout << "\nFinal result k:\n";
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < s; j++)
			cout << "  " << k[i][j];
		cout << endl;
	}
}

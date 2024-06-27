#ifndef SOFTMAX_H
#define SOFTMAX_H

//#include <iostream>
//using std::cout;
//using std::endl;
#include <cmath>
using std::sqrt;
using std::log;
using std::fabs;

#include "../../../../Eigen/Eigen/Dense"
using namespace Eigen;

namespace softmax {
	// computation result code
	enum ResultCode {
		NO_ERRORS = 0,
		UNKNOWN_PROBLEM,
		INCORRECT_SIZES,
		ICORRECT_HYPER_PARAMS,
	};

	template<typename Float, typename Int>
	int fitSoftmax(Float * features, Float * avgs, Float * stdevs, Int * target, 
		const int samplesCount, const int featuresCount, const int classesCount,
		int iterCount, Float learningRate, Float penalty, float tolerance,
		Float * params) noexcept {
		const int n = featuresCount;
		const int m = samplesCount;
		const int c = classesCount;

		Matrix<Float, Dynamic, Dynamic, ColMajor> X(n, m);

		Map<Matrix<Float, Dynamic, Dynamic, ColMajor>> XT(features, m, n);

		Map<Vector<Float, Dynamic>> avgX(avgs, n);
		Map<Vector<Float, Dynamic>> stdDevX(stdevs, n);

		Matrix<Float, Dynamic, Dynamic, ColMajor> Y(c, m);

		//Map<Matrix<Float, Dynamic, Dynamic, RowMajor>> W(params, c, n + 1);
		Matrix<Float, Dynamic, Dynamic, RowMajor> W(c, n);
		Vector<Float, Dynamic> B(c);

		Matrix<Float, Dynamic, Dynamic, ColMajor> Z(c, m);

		Matrix<Float, Dynamic, Dynamic, RowMajor> dZ(c, m);

		Matrix<Float, Dynamic, Dynamic, RowMajor> dW(c, n);
		Vector<Float, Dynamic> dB(c);

		Vector<Float, Dynamic> sums1(m);
		Vector<Float, Dynamic> sums2(c);

		uint32_t* classWeights = new uint32_t[c];
		Float loss = 0;
		Float lossPrev = 0;

		// 1. Prepare data

		// 1.1) Center features
		XT = XT.rowwise() - avgX.transpose();

		// 1.2) Scale features
		for (int i = 0; i < n; ++i)
			if (stdevs[i] > static_cast<Float>(0))
				XT.col(i) = XT.col(i) / stdDevX(i);

		// 1.3) Transpose
		X = XT.transpose();

		//cout << "X transposed:\n" << XT << endl;

		//cout << "\nX:\n" << X << endl;
		
		// 1.4) One-hot labels & class weights
		Y = Matrix<Float, Dynamic, Dynamic, ColMajor>::Zero(c, m);
		for (int i = 0; i < c; ++i)
			classWeights[i] = 0;

		for (int i = 0; i < m; ++i) {
			Y(target[i], i) = static_cast<Float>(1);
			classWeights[target[i]] += 1;
		}

        /*cout << "targets: ";
		for (int i = 0; i < m; ++i)
		  cout << target[i] << ", ";
		cout << endl;*/

		// 2. Init parameters		
		//W = Matrix<Float, Dynamic, Dynamic, RowMajor>::Random(c, n + 1);
		//W.col(n) = Vector<Float, Dynamic>::Zero(c);
		Float xavierScale = sqrt(6.0 / (c + n));
		W = Matrix<Float, Dynamic, Dynamic, RowMajor>::Random(c, n) * xavierScale;
		B = Vector<Float, Dynamic>::Zero(c);
		//cout << "\nW:\n" << W << endl;

		for (int iter = 0; iter < iterCount; ++iter) {
			//Z = ((W.block(0, 0, c, n) * X).colwise() + W.col(n)).exp();
			Z = ((W * X).colwise() + B);/// .exp();
			//cout << "\nZ:\n" << Z << endl;

			for (int j = 0; j < m; ++j)
				for (int i = 0; i < c; ++i)
					Z(i, j) = exp(Z(i, j)) * classWeights[i];
			
			sums1 = Z.colwise().sum();

			for (int i = 0; i < m; ++i)
				Z.col(i) = Z.col(i) / sums1(i);
						
			//cout << "\nZ:\n" << Z << endl;			

			//cout << "\nsums:\n" << sums << endl;				

            loss = 0;

			for (int i = 0; i < m; ++i)
				loss += -log(Z(target[i], i));

			loss /= m;

			if ((fabs(loss - lossPrev) < tolerance) || (loss != loss))
				break;
			
			lossPrev = loss;

			dZ = Z - Y;
			
			//cout << "\ndZ:\n" << dZ << endl;			

			dB = dZ.rowwise().mean();

			//cout << "\ndB:\n" << dB << endl;

			dW = dZ * XT / m;

			//cout << "\ndW:\n" << dW << endl;

			W = W * (static_cast<Float>(1) - learningRate * penalty / m) - dW * learningRate;

			//cout << "\nW:\n" << W << endl;			

			B = B - dB * learningRate;

			//cout << "\nB:\n" << B << endl;
		} // for iter

		//cout << "\nW:\n" << W << endl;
		//cout << "\nB:\n" << B << endl;		

		Map<Matrix<Float, Dynamic, Dynamic, RowMajor>> P(params, c, n + 1);
		P.block(0, 0, c, n) = W;
		P.col(n) = B;

		//std::cout << "\nP:\n" << P << std::endl;

		delete[] classWeights;

		return NO_ERRORS;
	}; // fitSoftmax
}; // softmax

#endif // !SOFTMAX_H
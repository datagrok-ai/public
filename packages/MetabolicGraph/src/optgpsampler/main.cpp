#include "sampler.cpp"
#include <iostream>
#include <fstream>
#include <chrono>

const double TOL = 1e-6;

const double DEFAULT_LOW_BOUND = -10;
const double DEFAULT_UPPER_BOUND = 10;

const unsigned SEED = 42;
const unsigned THINNING = 50;
const int NPROJ = 8;
const int WARM_STEPS = 2000;
const double FEAS_TOL = 1e-8;
const int FEAS_MAX_ITER = 4000;

double samplerRandomTest(const int metabolitesCount = 2, const int reactionsCount = 4, const int samplesCount = 1000) {
	Eigen::VectorXd beq(metabolitesCount);
	beq.setConstant(0.0);

	Eigen::VectorXd lb(reactionsCount);
	lb.setConstant(DEFAULT_LOW_BOUND);

	Eigen::VectorXd ub(reactionsCount);
	ub.setConstant(DEFAULT_UPPER_BOUND);

	Eigen::MatrixXd S(metabolitesCount, reactionsCount);
	//S << -2, 1, 0, 1, 0, 1, -2, 1;
	S.setRandom(); // [-1, 1]

	//std::cout << "S=\n" << S << "\n";

	LinearPolytope P;
	P.A = Eigen::MatrixXd();
	P.b = Eigen::VectorXd();
	P.Aeq = S; P.beq = beq;
	P.lb = lb; P.ub = ub;

	OptGPOptions opt;
	opt.seed = SEED;
	opt.thinning = THINNING;
	opt.nproj = NPROJ;
	opt.warmupSteps = WARM_STEPS;
	opt.feasTol = FEAS_TOL;
	opt.feasMaxIter = FEAS_MAX_ITER;

	OptGPSampler sampler(P, opt);

	auto start = std::chrono::high_resolution_clock::now();

	auto samples = sampler.sample(samplesCount);
	
	auto end = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> elapsed = end - start; // In milliseconds
	
	/*for (int i = 0; i < samplesCount; ++i)
		std::cout << i << "-th: " << samples[i].transpose() << "\n";*/

	double mad = 0;

	for (int i = 0; i < samplesCount; ++i) {
		mad = (S * samples[i]).norm();
		if (mad > TOL) {
			throw std::runtime_error("Invalid sample: " + std::to_string(i) + ", too big deviation: " + std::to_string(mad));
		}

		for (int k = 0; k < reactionsCount; ++k) {
			if ( (samples[i](k) > ub(k)) || (samples[i](k) < lb(k)) ) {
				throw std::runtime_error("Invalid sample: value is out of range");
			}
		}
	}	

	return elapsed.count();
}

/*
   NOTES:
   
     1) in sampler.cpp, the path to the Eigne lib is manually specified (see line 32)
   
     2) command for gcc: 
	     g++ -O3 -std=c++17 sampler.cpp main.cpp -o check-performance
*/

int main() {
	/*samplerRandomTest(50, 50, 1000);
	return 0;*/
	int reactMin = 0;
	int reactMax = 0;

	int metabMin = 0;
	int metabMax = 0;

	int step = 0;

	int samples = 0;

	std::cout << "Reactions min: ";
	std::cin >> reactMin;

	std::cout << "Reactions max: ";
	std::cin >> reactMax;

	std::cout << "Metabolites min: ";
	std::cin >> metabMin;

	std::cout << "Metabolites max: ";
	std::cin >> metabMax;

	std::cout << "Step: ";
	std::cin >> step;

	std::cout << "Samples: ";
	std::cin >> samples;

	std::ofstream file("eigen-" + std::to_string(samples) + "-samples-ms.csv", std::ios::out);

	file << "metabolites\\reactions";

	for (int i = reactMin; i <= reactMax; i += step)
		file << ',' << i;

	file << "\n";

	for (int metabolitesCount = metabMin; metabolitesCount <= metabMax; metabolitesCount += step) {
		file << metabolitesCount;

		for (int reactionsCount = reactMin; reactionsCount <= reactMax; reactionsCount += step) {
			std::cout << metabolitesCount << " x " << reactionsCount << ": ";
			file << ",";

			try {
				double execTime = samplerRandomTest(metabolitesCount, reactionsCount, samples);
				std::cout << execTime << " ms" << " --- " << "SUCCESS!\n";
				file << std::ceil(execTime);
			} catch (const std::runtime_error& e) {
				std::cerr << "Error: " << e.what() << std::endl;				
			}
			catch (...) { // Catch-all for any other unhandled exceptions
				std::cerr << "An unknown error occurred." << std::endl;
			}
		}

		file << "\n";
	}

	std::cout << "\nPress any key...";

	//
	//// Simple chain Sv=0 with bounds 0<=v<=10
	//Eigen::MatrixXd S(2, 4); S << -2, 1, 0, 1, 0, 1, -2, 1; // Aeq
	//Eigen::VectorXd beq(2); beq << 0, 0;
	//Eigen::VectorXd lb(4); lb << -10, -10, -10, -10;
	//Eigen::VectorXd ub(4); ub << 10, 10, 10, 10;

	//LinearPolytope P;
	//P.A = Eigen::MatrixXd();
	//P.b = Eigen::VectorXd();
	//P.Aeq = S; P.beq = beq;
	//P.lb = lb; P.ub = ub;

	//OptGPOptions opt; 
	//opt.seed = 42;
	//opt.thinning = 50;
	//opt.nproj = 8;
	//opt.warmupSteps = 2000;
	//opt.feasTol = 1e-8;
	//opt.feasMaxIter = 4000;

	//OptGPSampler sampler(P, opt);
	//auto samples = sampler.sample(1000);
	//std::cout << "N=" << samples.size() << " first sample: " << samples[0].transpose() << "\n";

	//std::cout << S * samples[0];

	std::cin.get();
}

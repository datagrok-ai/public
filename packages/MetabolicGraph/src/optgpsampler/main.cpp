#include "sampler.cpp"
#include <iostream>

int main() {
	// Simple chain Sv=0 with bounds 0<=v<=10
	Eigen::MatrixXd S(2, 3); S << 10, -2, 0, 0, 2, -10; // Aeq
	Eigen::VectorXd beq(2); beq << 0, 0;
	Eigen::VectorXd lb(3); lb << 0, 0, 0;
	Eigen::VectorXd ub(3); ub << 10, 10, 10;

	LinearPolytope P;
	P.A = Eigen::MatrixXd();
	P.b = Eigen::VectorXd();
	P.Aeq = S; P.beq = beq;
	P.lb = lb; P.ub = ub;

	OptGPOptions opt; opt.seed = 42; opt.thinning = 50; opt.nproj = 8; opt.warmupSteps = 2000;

	OptGPSampler sampler(P, opt);
	auto samples = sampler.sample(100);
	std::cout << "N=" << samples.size() << " first sample: " << samples[0].transpose() << "\n";

	std::cout << S * samples[0];

	std::cin.get();
}
